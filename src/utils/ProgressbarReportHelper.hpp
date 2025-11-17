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
        std::cout << " %\r";
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
        std::cout << "] " << static_cast<std::size_t>(progress * 100.0);
    }

    inline void progresstext(const std::string &str, const std::size_t textWidth = 20)
    {
        if (str.size() > textWidth)
            std::cout << str.substr(0, textWidth - 3) << "...";
        else
            std::cout << str;
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

    /**
     * A progressbar helper. API is simple. If subtasks are wanted, subtask() should be called before delegating.
     */
    struct Progressbar {
        Progressbar():m_reportlevel(ProgressbarReportLevel::FullReport){}
        Progressbar(ProgressbarReportLevel level):m_reportlevel(level){}

        void start(const std::string& str){
            ++m_processlevel;
            if(m_tasklevel != m_processlevel - 1){
                return;
            }

            // check if we already got a name
            if(m_names[m_tasklevel].empty()){
                m_names[m_tasklevel] = str;
            }
            if(m_tasklevel == 0){
                m_values[m_tasklevel].target = 1.0f;
                m_begin = std::chrono::steady_clock::now();
            }
            m_values[m_tasklevel].current = 0.0f;
            m_infos[m_tasklevel] = "";

            report();
        }
        void end(){
            if(m_processlevel == 0)
            {
                // we cannot end something which did not start
                return;
            }
            --m_processlevel;
            if(m_tasklevel != m_processlevel){
                return;
            }

            // whole task finished?
            if(m_tasklevel == 0){
                finish_report();
            }

            // delete old info and name
            m_names[m_tasklevel] = "";
            m_infos[m_tasklevel] = "";

            if(m_tasklevel == 0){
                return;
            }
            --m_tasklevel;
            m_values[m_tasklevel].current += m_values[m_tasklevel + 1].target;
            report();
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

        void update(float progress){
            if(m_tasklevel != m_processlevel - 1){
                return;
            }

            m_values[m_tasklevel].current = progress;
            report();
        }
        void update(const std::string& str){
            if(m_tasklevel != m_processlevel - 1){
                return;
            }

            m_infos[m_tasklevel] = str;
            report();
        }

        /**
         * Create a subtask.
         */
        void subtask(float completion){
            // no deep subtasks
            if(m_maxtasks > m_tasklevel + 1){
                return;
            }
            ++m_tasklevel;
            m_values[m_tasklevel].target = completion;
        }

        /**
         * Create a subtask with custom name
         */
        void subtask(const std::string& str, float completion){
            // no deep subtasks
            if(m_maxtasks > m_tasklevel + 1){
                return;
            }
            ++m_tasklevel;
            m_names[m_tasklevel] = str;
            m_values[m_tasklevel].target = completion;
        }

        // level does not need to be changed at runtime
        // void level(ProgressbarReportLevel reportLevel){
        //     m_reportlevel = reportLevel;
        // }
        
        private:
        /**
         * Reporting is done here (we delete all lines and recreate them)
         */
        void report(){
            // debounce (1s)
            if(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - m_begin).count() < 1000){
                return;
            }

            // update overall values (if necessary)
            if(m_reportlevel == ProgressbarReportLevel::FullReport)
            {
                m_values[m_tasklevel].overall = m_values[m_tasklevel].current;
                for(uint8_t i = 1; i <= m_tasklevel; ++i){
                    uint8_t j = m_tasklevel - i;
                    m_values[j].overall = m_values[j].current + m_values[j+1].target * m_values[j+1].overall;
                }
            }

            // report progress update
            erase_lines(m_drawnlines);

            switch (m_reportlevel)
            {
            case ProgressbarReportLevel::NoReport:
            case ProgressbarReportLevel::TimeReport:
                m_drawnlines = 0;
            case ProgressbarReportLevel::BusyReport:
                for(uint8_t i = 0; i <= m_tasklevel; ++i){
                    display(m_names[i], m_infos[i]);
                }
                m_drawnlines = m_tasklevel + 1;
                break;
            case ProgressbarReportLevel::FullReport:
            default:
                for(uint8_t i = 0; i <= m_tasklevel; ++i){
                    display(m_names[i], m_values[i].overall, m_infos[i]);
                }
                m_drawnlines = m_tasklevel + 1;
                break;
            }
        }

        void finish_report()
        {
            erase_lines(m_drawnlines);
            if(m_reportlevel != ProgressbarReportLevel::NoReport){
                display_finish(m_names[m_tasklevel], std::chrono::steady_clock::now() - m_begin);
            }
        }


        private:
        // How much information should be displayed
        ProgressbarReportLevel m_reportlevel;
        // current stack index
        std::uint8_t m_tasklevel = 0;
        // current level (count of start-end pairs)
        std::uint8_t m_processlevel = 0;
        // number of lines drawn
        std::uint8_t m_drawnlines = 0;
        // max (sub)tasks
        std::uint8_t m_maxtasks = 5;
        // stack of process/stage names
        std::array<std::string, 5> m_names;
        // further information for each process/stage
        std::array<std::string, 5> m_infos;
        // values for each process/stage
        std::array<ProgressValue, 5> m_values;
        // time when the main task was started
        std::chrono::steady_clock::time_point m_begin;
    };

    struct ProgressbarReportDynamic {
        ProgressbarReportDynamic():m_level(ProgressbarReportLevel::BusyReport){}
        ProgressbarReportDynamic(ProgressbarReportLevel level):m_level(level){}

        void start(const std::string& str){
            m_names.push_back(str);
            m_values.push_back(0.0f);
            m_infos.push_back("");
            if(m_level == ProgressbarReportLevel::NoReport){
                return;
            }

            if(m_names.size() == 1){
                m_begin = std::chrono::steady_clock::now();
            }
            else{
                // create new subtask
                std::cout << std::endl;
            }

            if (m_level == ProgressbarReportLevel::BusyReport){
                display(m_names.back(), m_infos.back());
            }
            else if (m_level == ProgressbarReportLevel::FullReport){
                display(m_names.back(), 0.0f, m_infos.back());
            }
        }
        void end(){
            if(m_names.size() == 0){
                return;
            }
            // main task?
            if(m_names.size() == 1){
                if (m_level != ProgressbarReportLevel::NoReport) {
                    erase_lines();
                    display_finish(m_names.back(), std::chrono::steady_clock::now() - m_begin);
                }
            }
            // pop subtask
            m_names.pop_back();
            m_values.pop_back();
            m_infos.pop_back();

            if(m_names.size() > 0){
                if (m_level == ProgressbarReportLevel::BusyReport){
                    erase_lines(2);
                    display(m_names.back(), m_infos.back());
                }
                else if(m_level == ProgressbarReportLevel::FullReport){
                    erase_lines(2);
                    display(m_names.back(), m_values.back(), m_infos.back());
                }
            }
        }
        bool stop(){
            return false;
        }
        void update(float progress){
            m_values.back() = progress;
            if (m_level == ProgressbarReportLevel::FullReport){
                erase_lines();
                display(m_names.back(), m_values.back(), m_infos.back());
            }
        }
        void update(const std::string& str){
            m_infos.back() = str;
            if (m_level == ProgressbarReportLevel::BusyReport){
                erase_lines();
                display(m_names.back(), m_infos.back());
            }
            else if(m_level == ProgressbarReportLevel::FullReport){
                erase_lines();
                display(m_names.back(), m_values.back(), m_infos.back());
            }
        }

        void level(ProgressbarReportLevel reportLevel){
            m_level = reportLevel;
        }

        private:
        ProgressbarReportLevel m_level;
        std::vector<std::string> m_names;
        std::vector<float> m_values;
        std::vector<std::string> m_infos;
        std::chrono::steady_clock::time_point m_begin;
    };
}
