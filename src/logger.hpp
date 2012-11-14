/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef ISOLATOR_LOGGER_HPP
#define ISOLATOR_LOGGER_HPP


#include <boost/thread.hpp>
#include <boost/timer/timer.hpp>
#include <cstdio>
#include <deque>
#include <queue>
#include <string>


/* This file consists of some overly-elaborate code giving feedback on the
 * terminal, including progress bars and logging. */


/* This class contains necessary details for a task, which the logger will
 * report the progress of.
 *
 * LoggerTasks form a tree: a LoggerTask may have subtasks.
 *
 * This class is not constructed directly, but rather by calling
 * `Logger::push_task`.
 */
class LoggerTask
{
    public:
        LoggerTask(const LoggerTask&);
        ~LoggerTask();

        /* Add a subtask.
         *
         * Args:
         *   name: Name of the sub-task.
         *   n: Number of work units.
         *   max_subtask_show: At most this many subtasks are shown.
         * */
        void push_task(const char* name, size_t n = 0, size_t max_subtask_show = 0);

        /* Remove the subtask whith the given name. */
        void pop_task(const char* name);

        /* Get the subtask with the given name. */
        LoggerTask& get_task(const char* name);

        /* Update with the completion of one unit of work. */
        void inc(size_t d = 1);

    private:
        boost::timer::cpu_timer timer;

        void unset_drawn();
        void print(int indent = 0);
        bool clear();

        LoggerTask();
        LoggerTask(const char* name, size_t n = 0, size_t max_subtask_show = 0);

        std::string name;
        size_t k, n;

        size_t max_subtask_show;

        /** Has the task beet printed yet. */
        bool drawn;

        /** Tasks must be poped/deleted lazily, since they are responsible for
         * clearing what they have written.  This is set when the task has been
         * popped from its parent task. */
        bool popped;

        /* Child sub-tasks indexed by name. */
        std::map<std::string, LoggerTask> children;

        boost::mutex mut;

        friend class Logger;
};


/**
 * A singleton log writing class.
 *
 */
class Logger
{
    public:
        static Logger& instance();
        static void start();
        static void end();

        ~Logger();

        /* Logger levels. */
        enum level
        {
            DEBUG,
            INFO,
            WARN,
            ERROR
        };

        /* Output log info at the given level. */
        static void debug (const char* fmt, ...);
        static void info  (const char* fmt, ...);
        static void warn  (const char* fmt, ...);
        static void error (const char* fmt, ...);

        /* print a message and exit. */
        static void abort (const char* fmt, ...);
        static void abort (const char* fmt, va_list);

        /* Print a log message at the given level. */
        static void put(level, const char* fmt, ...);

        /* More low-level printing. No filtering, etc. */
        static void print(const char* msg);

        /* Set the logger level. Messages at or above the logger level will be
         * printed. */
        static void set_level(level);

        /* Add a task. */
        static void push_task(const char* name, size_t num_steps = 0,
                              size_t max_subtask_show = 0);

        /* Remove a task by name. */
        static void pop_task(const char* name);

        /* Get a task by name. */
        static LoggerTask& get_task(const char* name);

        /* Stop printing output until `reume` is called. */
        static void suspend();

        /* Resume printing output after a call to suspend. */
        static void resume();

        /* Print any buffered output. */
        void flush();

        /* True if output should be in color. */
        bool color;

    private:
        Logger();
        Logger(Logger const&);
        void operator = (Logger const&);

        /* Logger is suspended. No printing. */
        bool suspended;

        /* Print a formatted string at the given level. */
        void put(level, const char* fmt, va_list);

        /* Output file. */
        FILE* fout;

        /* Logger level. */
        level L;

        /* True when logger is finished and should terminate as soon as
         * possible. */
        bool finished;

        /* Run the print_loop function in the background. */
        boost::thread* print_thread;

        /* An infinite loop that runs in the background printing progress bars.
         * and other logger output. */
        void print_loop();

        /* Messages to be printed. */
        std::queue<std::pair<level, std::string> > log_queue;

        /* The root task, to which all other are children. */
        LoggerTask root_task;

        /* Some space used for string formatting. */
        char* msg_buffer;
        static const size_t msg_buffer_len;

        boost::mutex mut;
};



#endif



