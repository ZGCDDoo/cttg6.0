#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "MPIUtilities.hpp"

namespace Logging
{

const std::string ROOT = "ROOT";

void Init(const std::string &loggerName = ROOT, const std::string logLevel = "DEBUG")
{
    // #ifdef HAVEMPI
    // mpi::environment env;
    // mpi::communicator world;
    // #endif

    if (mpiUt::Rank() == mpiUt::master)
    {
        auto logger_ = spdlog::stdout_color_st(loggerName);

        if (logLevel == "INFO")
        {
            logger_->set_level(spdlog::level::info);
        }
        else if (logLevel == "WARN")
        {
            logger_->set_level(spdlog::level::warn);
        }
        else if (logLevel == "DEBUG")
        {
            logger_->set_level(spdlog::level::debug);
        }
        else
        {
            logger_->set_level(spdlog::level::critical);
        }

        const std::string msg = "Logger " + loggerName + " initialized";
        logger_->info(msg);
    }
}

// template <typename... TArgs>
void Info(const std::string &msg, const std::string loggerName = ROOT)
{
    // #ifdef HAVEMPI
    // mpi::environment env;
    // mpi::communicator world;
    // #endif
    if (mpiUt::Rank() == mpiUt::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            return; //throw std::runtime_error("Ayaya, logger not found. Stupido !");
        }
        logger->info(msg);
    }
}

void Warn(const std::string &msg, const std::string loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Rank() == mpiUt::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            return; //throw std::runtime_error("Ayaya, logger not found. Stupido !");
        }

        logger->warn(msg);
    }
}

void Debug(const std::string &msg, const std::string loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Rank() == mpiUt::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            return; //throw std::runtime_error("Ayaya, logger not found. Stupido !");
        }
        logger->debug(msg);
    }
}

void Critical(const std::string &msg, const std::string loggerName = ROOT)
{
    // #ifdef HAVEMPI
    //     mpi::environment env;
    //     mpi::communicator world;
    // #endif
    if (mpiUt::Rank() == mpiUt::master)
    {
        auto logger = spdlog::get(loggerName);
        if (!logger)
        {
            return; //throw std::runtime_error("Ayaya, logger not found. Stupido !");
        }
        logger->critical(msg);
    }
}

} // namespace Logging
