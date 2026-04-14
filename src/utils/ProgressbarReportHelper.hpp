#pragma once

#include <string>
#include <iostream>
#include <chrono>
#include <vector>
#include <array>

namespace progressbar
{
    inline void flush()
    {
        std::cout << "\r";
        std::cout.flush();
    }

    inline void erase_lines(std::size_t count = 1)
    {
        if (count == 0)
        {
            return;
        }
        std::cout << "\x1b[2K"; // Delete current line
        for (std::size_t i = 1; i < count; i++)
        {
            std::cout
                << "\x1b[1A"  // Move cursor up one
                << "\x1b[2K"; // Delete the entire line
        }
        flush();
    }

    inline void progressbar(float progress, const std::size_t barWidth = 70)
    {
        std::cout << "[";
        std::size_t pos = static_cast<std::size_t>(barWidth * progress);
        for (std::size_t i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] " << static_cast<std::size_t>(progress * 100.0) << " %";
    }

    inline void progressbar(std::stringstream& output, float progress, const std::size_t barWidth = 70)
    {
        output << "[";
        std::size_t pos = static_cast<std::size_t>(barWidth * progress);
        for (std::size_t i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                output << "=";
            else if (i == pos)
                output << ">";
            else
                output << " ";
        }
        output << "] " << static_cast<std::size_t>(progress * 100.0) << " %";
    }

    inline void progresstext(const std::string &text, const std::size_t textWidth = 20)
    {
        if (text.size() > textWidth)
            std::cout << text.substr(0, textWidth - 3) << "...";
        else
            std::cout << text;
    }

    inline void progresstext(std::stringstream& output, const std::string &text, const std::size_t textWidth = 20)
    {
        if (text.size() > textWidth)
            output << text.substr(0, textWidth - 3) << "...";
        else
            output << text;
    }

    inline void display(const std::string &progressname, float progress, const std::string &stepname = "")
    {
        progresstext(progressname, 20);
        std::cout << "   " << stepname;
        progressbar(progress, 70);
        flush();
    }

    inline void display(const std::string &progressname, const std::string &stepname)
    {
        progresstext(progressname, 20);
        std::cout << "   " << stepname;
        flush();
    }

    inline void display_finish(const std::string &progressname)
    {
        progresstext(progressname, 90);
        std::cout << " finished!" << std::endl;
    }

    template <class D>
    inline void display_finish(const std::string &progressname, const D&time){
        progresstext(progressname, 90);
        std::cout << " finished! [" << (std::chrono::duration_cast<std::chrono::microseconds>(time).count() / 1000000.0) << " sec]" << std::endl;
    }

    enum ProgressbarReportLevel {
        // Do not report anything
        NoReport,
        // Only report finished process
        TimeReport,
        // Report current step
        BusyReport,
        // Report current step and progress
        FullReport
    };

    struct ProgressValue {
        // how big the target value is (conversion factor to parent task)
        float target;
        // current value (ignoring subtasks)
        float current;
        // buffer to read (overall progress, including subtasks)
        float overall;
    };

    struct TaskLevelData {
        std::chrono::steady_clock::time_point begin;
        std::string name;
        std::string info;

        void display()
        {
            progresstext(name, 90 - info.size());
            std::cout << "   " << info;
            flush();
        }
    };

    /**
     * A progressbar helper. API is simple.
     */
    struct Reporter {
        Reporter():m_reportlevel(ProgressbarReportLevel::FullReport){}
        Reporter(ProgressbarReportLevel level):m_reportlevel(level){}

        ProgressbarReportLevel level(){
            return m_reportlevel;
        }

        void level(ProgressbarReportLevel level){
            m_reportlevel = level;
        }

        void start(const std::string& str){
            ++m_tasklevel;
            // // start main clock if not done yet
            // if(m_tasklevel == 0 && std::chrono::duration_cast<std::chrono::milliseconds>(m_begin.time_since_epoch()).count() == 0)
            // {
            //     m_begin = std::chrono::steady_clock::now();
            // }
            if(m_tasklevel > 0){
                return;
            }
            m_task.begin = std::chrono::steady_clock::now();
            m_task.name = str;
            m_task.info = "";
            report();
        }
        void end(){
            // task finished?
            if(m_tasklevel == 0){
                finish_report();
                // delete old info and name
                m_task.name = "";
                m_task.info = "";   
            }
            
            --m_tasklevel;

        }
        bool stop(){
            return false;
        }
        bool stop_and_end(){
            bool stopped = stop();
            if(stopped){
                end();
            }
            return stopped;
        }

        void update(const std::string& info){
            if(m_tasklevel > 0){
                return;
            }

            m_task.info = info;
            report();
        }

        void update(const std::string& info, const float progress){
            if(m_tasklevel > 0){
                return;
            }

            std::stringstream output;
            output << info;
            progressbar(output, progress);

            m_task.info = output.str();
            report();
        }

        void update(const float progress){
            if(m_tasklevel > 0){
                return;
            }

            std::stringstream output;
            progressbar(output, progress);

            m_task.info = output.str();
            report();
        }
        
        private:
        /**
         * Reporting is done here (we delete all lines and recreate them)
         */
        void report(){
            // report progress update
            // erase_lines(1);

            switch (m_reportlevel)
            {
            case ProgressbarReportLevel::NoReport:
            case ProgressbarReportLevel::TimeReport:
            case ProgressbarReportLevel::BusyReport:
            case ProgressbarReportLevel::FullReport:
            default:
                m_task.display();
                break;
            }
        }

        void finish_report()
        {
            // erase_lines(1);
            if(m_reportlevel != ProgressbarReportLevel::NoReport){
                display_finish(m_task.name, std::chrono::steady_clock::now() - m_task.begin);
            }
        }


        private:
        // How much information should be displayed
        ProgressbarReportLevel m_reportlevel;
        // current task depth
        std::uint8_t m_tasklevel = 255;
        TaskLevelData m_task;
    };
} // namespace progressbar
