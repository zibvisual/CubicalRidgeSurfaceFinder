#pragma once

#include <string>
#include <iostream>
#include <chrono>
#include <vector>

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

    // clang-format off
    struct ProgressbarReportNone {
        void start(const std::string& str){/* noop */}
        void end(){/* noop */}
        bool stop(){return false;}
        void update(float){/* noop */}
        void update(const std::string& str){/* noop */}
    };

    struct ProgressbarReportBusy {
        void start(const std::string& str){
            m_progress = str;
            display(m_progress, "");
        }
        void end(){
            erase_lines();
            display_finish(m_progress);
        }
        bool stop(){return false;}
        void update(float){/* noop */}
        void update(const std::string& str){
            erase_lines();
            display(m_progress, str);
        }

        private:
        std::string m_progress;
    };

    struct ProgressbarReportFull {
        void start(const std::string& str){
            m_progress = str;
            m_info = "";
            m_value = 0.0f;
            display(m_progress, m_value, m_info);
            m_begin = std::chrono::steady_clock::now();
        }
        void end(){
            erase_lines();
            display_finish(m_progress, std::chrono::steady_clock::now() - m_begin);
        }
        bool stop(){return false;}
        void update(float progress){
            m_value = progress;
            erase_lines();
            display(m_progress, m_value, m_info);
        }
        void update(const std::string& str){
            m_info = str;
            erase_lines();
            display(m_progress, m_value, m_info);
        }

        private:
        std::string m_progress;
        std::string m_info;
        float m_value;
        std::chrono::steady_clock::time_point m_begin;
    };

    // clang-format on

    enum ProgressbarReportLevel {
        NoReport,
        BusyReport,
        FullReport
    };

    struct SimpleProgressbarReportDynamic {
        SimpleProgressbarReportDynamic():reportLevel(ProgressbarReportLevel::BusyReport){}
        SimpleProgressbarReportDynamic(ProgressbarReportLevel level):reportLevel(level){}

        void start(const std::string& str){
            m_progress = str;
            m_info = "";
            m_value = 0.0f;
            if (reportLevel == ProgressbarReportLevel::BusyReport) display(m_progress, m_value, m_info);
            else if (reportLevel == ProgressbarReportLevel::FullReport) display(m_progress, m_value, m_info);
            if(timer){
                m_begin = std::chrono::steady_clock::now();
            }
        }
        void end(){
            if (reportLevel != ProgressbarReportLevel::NoReport) {
                erase_lines();
                if (timer){
                    display_finish(m_progress, std::chrono::steady_clock::now() - m_begin);
                }else{
                    display_finish(m_progress);
                }
            }
        }
        bool stop(){
            return false;
        }
        void update(float progress){
            m_value = progress;
            if (reportLevel == ProgressbarReportLevel::FullReport){
                erase_lines();
                display(m_progress, m_value, m_info);
            }
        }
        void update(const std::string& str){
            m_info = str;
            if (reportLevel == ProgressbarReportLevel::BusyReport){
                erase_lines();
                display(m_progress, m_value, m_info);
            }
        }

        ProgressbarReportLevel reportLevel;
        bool timer = true;

        private:
        std::string m_progress;
        std::string m_info;
        float m_value;
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
