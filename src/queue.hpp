
#ifndef ISOLATOR_QUEUE_HPP
#define ISOLATOR_QUEUE_HPP

#include <queue>
#include <boost/thread.hpp>

/*
 * A producer/consumer queue for parallel processing.
 */
template <typename T>
class Queue
{
    public:
        Queue(size_t max_size = 0)
            : max_size(max_size)
        {}

        ~Queue() {}

        /*
         * Push a value onto the queue.
         */
        void push(const T& x)
        {
            {
                boost::unique_lock<boost::mutex> lock(mut);

                while (max_size > 0 && q.size() >= max_size) {
                    pop_cond.wait(lock);
                }

                q.push(x);
            }

            push_cond.notify_one();
        }


        /*
         * Pop a value from the queue. If there are none, the thread sleeps
         * until one is available.
         */
        T pop()
        {
            T v;

            {
                boost::unique_lock<boost::mutex> lock(mut);

                while (q.empty()) {
                    push_cond.wait(lock);
                }

                v = q.front();
                q.pop();
            }

            pop_cond.notify_one();
            return v;
        }


        bool empty()
        {
            return q.empty();
        }


    private:
        size_t max_size;
        boost::condition_variable push_cond;
        boost::condition_variable pop_cond;
        boost::mutex mut;
        std::queue<T> q;
};

#endif

