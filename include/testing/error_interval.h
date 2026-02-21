#ifndef ERROR_INTERVAL_H
#define ERROR_INTERVAL_H

struct ErrorInterval {
    const double m_value;
    const double m_error;
    ErrorInterval(double value, double error) : m_value(value), m_error(error) {}
};

#endif
