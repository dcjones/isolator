
#include <cstdarg>
#include <cassert>
#include <ctime>
#include <cstdio>

#include "config.h"
#include "logger.hpp"


Logger& Logger::instance()
{
    static Logger S;
    return S;
}

const size_t Logger::msg_buffer_len = 1024;


Logger::Logger()
    : color(true)
    , suspended(false)
    , fout(stderr)
    , L(INFO)
    , print_thread(NULL)
    , msg_buffer(new char[msg_buffer_len])
{
}


Logger::~Logger()
{
    finished = true;
    if (print_thread) {
        delete print_thread;
    }
    fflush(stdout);
    delete [] msg_buffer;
}


void Logger::start()
{
    boost::lock_guard<boost::mutex> lock(instance().mut);

    if (instance().print_thread == NULL) {
        instance().print_thread =
            new boost::thread(boost::bind(&Logger::print_loop, &instance()));
    }
}


void Logger::end()
{
    instance().flush();
    boost::lock_guard<boost::mutex> lock(instance().mut);

    instance().finished = true;
    instance().print_thread->join();

    delete instance().print_thread;
    instance().print_thread = NULL;
}


void Logger::set_level(level L)
{
    instance().L = L;
}


void Logger::debug(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    Logger::instance().put(DEBUG, fmt, args);

    va_end(args);
}


void Logger::info(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    Logger::instance().put(INFO, fmt, args);

    va_end(args);
}


void Logger::warn(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    Logger::instance().put(INFO, fmt, args);

    va_end(args);
}


void Logger::error(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    Logger::instance().put(ERROR, fmt, args);

    va_end(args);
}


void Logger::abort(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    Logger::instance().put(ERROR, fmt, args);
    va_end(args);

    if (Logger::instance().suspended) resume();
    Logger::instance().flush();

    exit(EXIT_FAILURE);
}


void Logger::abort(const char* fmt, va_list args)
{
    Logger::instance().put(INFO, fmt, args);
    exit(EXIT_FAILURE);
}


void Logger::put(level L, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    Logger::instance().put(L, fmt, args);

    va_end(args);
}


void Logger::put(level l, const char* fmt, va_list args)
{
    boost::lock_guard<boost::mutex> lock(mut);

    int len = vsnprintf(msg_buffer, msg_buffer_len - 1, fmt, args);

    /* make sure there are no trailing newlines */
    while (len > 0 && msg_buffer[len - 1] == '\n') {
        msg_buffer[len - 1] = '\0';
        --len;
    }

    log_queue.push(std::pair<level, std::string>(l, msg_buffer));
}


void Logger::print(const char* msg)
{
    fputs(msg, stdout);
    fflush(stdout);
}


void Logger::print_loop()
{
    while (!finished)
    {
        if (!suspended) flush();
        boost::this_thread::sleep(boost::posix_time::milliseconds(250));
    }
}


void Logger::suspend()
{
    Logger::instance().flush();
    Logger::instance().suspended = true;
}


void Logger::resume()
{
    Logger::instance().root_task.unset_drawn();
    Logger::instance().suspended = false;
}


void Logger::flush()
{
    boost::lock_guard<boost::mutex> lock(mut);

    if (suspended) return;

    static time_t t;
    static struct tm ts;
    static char time_str[200];

    root_task.clear();

    if (!log_queue.empty()) {
        t = time(NULL);
        localtime_r(&t, &ts);
        // RFC-2822 date
        strftime(time_str, 200, "%a, %d %b %Y %T %z", &ts);
    }

    while (!log_queue.empty()) {
        if (log_queue.front().first >= L) {
            printf("[%s] ", time_str);

            if (color) {
                switch (log_queue.front().first) {
                    case DEBUG: printf("\033[37m"); break;
                    case INFO:  break;
                    case WARN:  printf("\033[33m"); break;
                    case ERROR: printf("\033[31m"); break;
                }
            }

            printf("%s\n", log_queue.front().second.c_str()) ;

            if (color) printf("\033[00m");
        }

        log_queue.pop();
    }

    root_task.print();
}


void Logger::push_task(const char* name, size_t num_steps, size_t max_subtask_show)
{
    boost::lock_guard<boost::mutex> lock(Logger::instance().mut);
    Logger::instance().root_task.push_task(name, num_steps, max_subtask_show);
}


void Logger::pop_task(const char* name)
{
    boost::lock_guard<boost::mutex> lock(Logger::instance().mut);
    Logger::instance().root_task.pop_task(name);
}


LoggerTask& Logger::get_task(const char* name)
{
    return Logger::instance().root_task.get_task(name);
}


LoggerTask::LoggerTask()
    : name()
    , k(0)
    , n(0)
    , max_subtask_show(0)
    , drawn(false)
    , popped(false)
{
}


LoggerTask::LoggerTask(const char* name, size_t n, size_t max_subtask_show)
    : name(name)
    , k(0)
    , n(n)
    , max_subtask_show(max_subtask_show)
    , drawn(false)
    , popped(false)
{
}


LoggerTask::LoggerTask(const LoggerTask& other)
    : name(other.name)
    , k(other.k)
    , n(other.n)
    , max_subtask_show(other.max_subtask_show)
    , drawn(false)
    , popped(false)
{
}


LoggerTask::~LoggerTask()
{
}


void LoggerTask::push_task(const char* name, size_t n, size_t max_subtask_show)
{
    boost::lock_guard<boost::mutex> lock(mut);

    children.insert(
        std::pair<std::string, LoggerTask>(
            name, LoggerTask(name, n, max_subtask_show)));
}


void LoggerTask::pop_task(const char* name)
{
    boost::lock_guard<boost::mutex> lock(mut);

    std::map<std::string, LoggerTask>::iterator i = children.find(name);

    if (i != children.end()) {
        i->second.popped = true;
    }
}


LoggerTask& LoggerTask::get_task(const char* name)
{
    std::map<std::string, LoggerTask>::iterator i = children.find(name);

    assert(i != children.end());

    return i->second;
}


void LoggerTask::inc(size_t d)
{
    boost::lock_guard<boost::mutex> lock(mut);

    if (k == 0) {
        timer.start();
    }
    if (k + d <= n || n == 0) k += d;
}


void LoggerTask::unset_drawn()
{
    boost::lock_guard<boost::mutex> lock(mut);
    drawn = false;
    for (std::map<std::string, LoggerTask>::reverse_iterator i = children.rbegin();
         i != children.rend();
         ++i) {
        i->second.unset_drawn();
    }
}


bool LoggerTask::clear()
{
    boost::lock_guard<boost::mutex> lock(mut);

    if (!drawn) return false;

    /* clear children */
    size_t j = 0;
    for (std::map<std::string, LoggerTask>::reverse_iterator i = children.rbegin();
         i != children.rend();
         ++i) {
        if (i->second.clear()) ++j;
    }

    printf("\033[A\033[2K\033[G"); /* move up, erase line, move to start */
    if (!name.empty() && n > 0) {
        printf("\033[A\033[2K\033[G");
        printf("\033[A\033[2K\033[G");
    }

    /* delete marked children, since we need them no longer */
    std::map<std::string, LoggerTask>::iterator i = children.begin();
    while (i != children.end()) {
        if (i->second.popped) {
            children.erase(i++);
        }
        else ++i;
    }

    return true;
}


void LoggerTask::print(int indent)
{
    boost::lock_guard<boost::mutex> lock(mut);

    if (!name.empty() && n == 0) {
        for (int i = 0; i < indent; ++i) printf("  ");
        printf("%s ...", name.c_str());
        for (size_t i = 0; i < k; ++i) {
            putchar('.');
        }
    }

    /* draw a fancy progress bar, if we know how many steps */
    else if (!name.empty()) {
        for (int i = 0; i < indent; ++i) printf("  ");

        printf("%s:\n", name.c_str());
        for (int i = 0; i < indent; ++i) printf("  ");
        printf("    ");

        if (Logger::instance().color) printf("\033[32m");

        printf("[");

        double step = 80.0 / (double) n;
        size_t i = 0;
        if ((int) ((double) k * step ) > 0) {
            for (; i < (size_t) ((double) k * step) - 1; ++i) {
                putchar('=');
            }
            putchar('>');
        }

        for (; i < 79; ++i) putchar(' ');

        printf("] ");

        if (Logger::instance().color) printf("\033[01m");

        /* percent completed */
        printf("%4.1f%%", 100.0 * (double) k / (double) n);

        /* time remaining */
        if (k <= 1) {
            printf("    ?:?? ETA\n");
        }
        else {
            boost::timer::cpu_times t = timer.elapsed();
            boost::timer::nanosecond_type d = t.wall / (k - 1);
            boost::timer::nanosecond_type remaining = (n - k) * d;
            double rem_sec = (double) remaining / 1e9;
            printf("%5lu:%02lu ETA\n",
                   (unsigned long) round(rem_sec / 60.0),
                   (unsigned long) round(fmod(rem_sec, 60.0)));
        }

        if (Logger::instance().color) printf("\033[00m");
    }
    putchar('\n');

    /* draw children */
    size_t j = 0;
    std::map<std::string, LoggerTask>::iterator i;
    for (i = children.begin();
         i != children.end();
         ++j, ++i) {

        if (i->second.popped) continue;

        i->second.print(indent + 1);
        ++j;
    }

    drawn = true;
}

