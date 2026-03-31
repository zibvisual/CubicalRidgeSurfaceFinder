#include "common.hpp"
  
// Definition FIX
const stype_t stype_t_helper<int8_t>::val;
const stype_t stype_t_helper<int16_t>::val;
const stype_t stype_t_helper<int32_t>::val;
const stype_t stype_t_helper<int64_t>::val;
const stype_t stype_t_helper<uint8_t>::val;
const stype_t stype_t_helper<uint16_t>::val;
const stype_t stype_t_helper<uint32_t>::val;
const stype_t stype_t_helper<uint64_t>::val;
const stype_t stype_t_helper<float>::val;
const stype_t stype_t_helper<double>::val;


// Parse Token Specialization
template <>
void
parseToken<unsigned char>(char* token, unsigned char& val)
{
    // cast to int first, won't work otherwise (reason unclear)
    try
    {
        int temp = std::stoi(token);
        val = static_cast<unsigned char>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<unsigned short>(char* token, unsigned short& val)
{
    try
    {
        int temp = std::stoi(token);
        val = static_cast<unsigned short>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<unsigned int>(char* token, unsigned int& val)
{
    try
    {
        int temp = std::stoi(token);
        val = static_cast<unsigned int>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<unsigned long long>(char* token, unsigned long long& val)
{
    try
    {
        int temp = std::stoi(token);
        val = static_cast<unsigned long long>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<char>(char* token, char& val)
{
    try
    {
        int temp = std::stoi(token);
        val = static_cast<char>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<short>(char* token, short& val)
{
    try
    {
        int temp = std::stoi(token);
        val = static_cast<short>(temp);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<int>(char* token, int& val)
{
    try
    {
        val = std::stoi(token);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<long long>(char* token, long long& val)
{
    try
    {
        val = std::stoll(token);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<float>(char* token, float& val)
{
    try
    {
        char* end;
        val = std::strtof(token, &end);
        if (end == token)
        {
            throw std::invalid_argument("Invalid argument");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}

template <>
void
parseToken<double>(char* token, double& val)
{
    try
    {
        char* end;
        val = std::strtod(token, &end);
        if (end == token)
        {
            throw std::invalid_argument("Invalid argument");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: Exception occurred - " << e.what() << std::endl;
    }
}