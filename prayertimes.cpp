#include "prayertimes.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>

double latitude;
double longitude;
double time_zone;
double julian_date;

MethodConfig method_params[CalculationMethodsCount];
CalculationMethod calc_method;        // caculation method
JuristicMethod asr_juristic;        // Juristic method for Asr
AdjustingMethod adjust_high_lats;    // adjusting method for higher latitudes
double dhuhr_minutes;        // minutes after mid-day for Dhuhr

/* ---------------------- Trigonometric Functions ----------------------- */

/* degree sin */
static double dsin(double d)
{
    return sin(deg2rad(d));
}

/* degree cos */
static double dcos(double d)
{
    return cos(deg2rad(d));
}

/* degree tan */
static double dtan(double d)
{
    return tan(deg2rad(d));
}

/* degree arcsin */
static double darcsin(double x)
{
    return rad2deg(asin(x));
}

/* degree arccos */
static double darccos(double x)
{
    return rad2deg(acos(x));
}

/* degree arctan */
static double darctan(double x)
{
    return rad2deg(atan(x));
}

/* degree arctan2 */
static double darctan2(double y, double x)
{
    return rad2deg(atan2(y, x));
}

/* degree arccot */
static double darccot(double x)
{
    return rad2deg(atan(1.0 / x));
}

/* degree to radian */
static double deg2rad(double d)
{
    return d * M_PI / 180.0;
}

/* radian to degree */
static double rad2deg(double r)
{
    return r * 180.0 / M_PI;
}

/* range reduce angle in degrees. */
static double fix_angle(double a)
{
    a = a - 360.0 * floor(a / 360.0);
    a = a < 0.0 ? a + 360.0 : a;
    return a;
}

/* range reduce hours to 0..23 */
static double fix_hour(double a)
{
    a = a - 24.0 * floor(a / 24.0);
    a = a < 0.0 ? a + 24.0 : a;
    return a;
}

/* ------------------------------------------------------------------------------------------------- */
/* calculate julian date from a calendar date */
double get_julian_date(int year, int month, int day)
{
    if (month <= 2)
    {
        year -= 1;
        month += 12;
    }

    double a = floor(year / 100.0);
    double b = 2 - a + floor(a / 4.0);

    return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + b - 1524.5;
}
/* compute prayer times at given julian date */
void compute_day_times(double times[])
{
    double default_times[] = { 5, 5, 5, 5, 5, 5, 5 };        // default times
    for (int i = 0; i < TimesCount; ++i)
        times[i] = default_times[i];

    for (int i = 0; i < NUM_ITERATIONS; ++i)
        compute_times(times);

    adjust_times(times);
}        

/* compute declination angle of sun and equation of time */
DoublePair sun_position(double jd)
{
    double d = jd - 2451545.0;
    double g = fix_angle(357.529 + 0.98560028 * d);
    double q = fix_angle(280.459 + 0.98564736 * d);
    double l = fix_angle(q + 1.915 * dsin(g) + 0.020 * dsin(2 * g));

    // double r = 1.00014 - 0.01671 * dcos(g) - 0.00014 * dcos(2 * g);
    double e = 23.439 - 0.00000036 * d;

    double dd = darcsin(dsin(e) * dsin(l));
    double ra = darctan2(dcos(e) * dsin(l), dcos(l)) / 15.0;
    ra = fix_hour(ra);
    double eq_t = q / 15.0 - ra;
    DoublePair dp = { dd, eq_t };
    return dp;
}

/* set the juristic method for Asr */
void set_asr_method(JuristicMethod method_id)
{
    asr_juristic = method_id;
}

/* set adjusting method for higher latitudes */
void set_high_lats_adjust_method(AdjustingMethod method_id)
{
    adjust_high_lats = method_id;
}

/* compute declination angle of sun */
double sun_declination(double jd)
{
    return sun_position(jd).first;
}

/* compute equation of time */
double equation_of_time(double jd)
{
    return sun_position(jd).second;
}

/* compute mid-day (Dhuhr, Zawal) time */
double compute_mid_day(double _t)
{
    double t = equation_of_time(julian_date + _t);
    double z = fix_hour(12 - t);
    return z;
}

/* compute time for a given angle G */
double compute_time(double g, double t)
{
    double d = sun_declination(julian_date + t);
    double z = compute_mid_day(t);
    double v = 1.0 / 15.0 * darccos((-dsin(g) - dsin(d) * dsin(latitude)) / (dcos(d) * dcos(latitude)));
    return z + (g > 90.0 ? -v : v);
}
/* compute prayer times at given julian date */
/* compute the time of Asr */
double compute_asr(int step, double t)  // Shafii: step=1, Hanafi: step=2
{
    double d = sun_declination(julian_date + t);
    double g = -darccot(step + dtan(fabs(latitude - d)));
    return compute_time(g, t);
}

void compute_times(double times[])
{
    day_portion(times);

    times[Fajr]    = compute_time(180.0 - method_params[calc_method].fajr_angle, times[Fajr]);
    times[Sunrise] = compute_time(180.0 - 0.833, times[Sunrise]);
    times[Dhuhr]   = compute_mid_day(times[Dhuhr]);
    times[Asr]     = compute_asr(1 + asr_juristic, times[Asr]);
    times[Sunset]  = compute_time(0.833, times[Sunset]);
    times[Maghrib] = compute_time(method_params[calc_method].maghrib_value, times[Maghrib]);
    times[Isha]    = compute_time(method_params[calc_method].isha_value, times[Isha]);
}
/* compute the difference between two times  */
static double time_diff(double time1, double time2)
{
    return fix_hour(time2 - time1);
}

static std::string int_to_string(int num)
{
    return std::to_string(num);
}

/* add a leading 0 if necessary */
static std::string two_digits_format(int num)
{
    char tmp[16];
    tmp[0] = '\0';
    sprintf(tmp, "%2.2d", num);
    return std::string(tmp);
}
/* the night portion used for adjusting times in higher latitudes */
double night_portion(double angle)
{
    switch (adjust_high_lats)
    {
        case AngleBased:
            return angle / 60.0;
        case MidNight:
            return 1.0 / 2.0;
        case OneSeventh:
            return 1.0 / 7.0;
        default:
            // Just to return something!
            // In original library nothing was returned
            // Maybe I should throw an exception
            // It must be impossible to reach here
            return 0;
    }
}
/* adjust Fajr, Isha and Maghrib for locations in higher latitudes */
void adjust_high_lat_times(double times[])
{
    double night_time = time_diff(times[Sunset], times[Sunrise]);        // sunset to sunrise

    // Adjust Fajr
    double fajr_diff = night_portion(method_params[calc_method].fajr_angle) * night_time;
    if (isnan(times[Fajr]) || time_diff(times[Fajr], times[Sunrise]) > fajr_diff)
        times[Fajr] = times[Sunrise] - fajr_diff;

    // Adjust Isha
    double isha_angle = method_params[calc_method].isha_is_minutes ? 18.0 : method_params[calc_method].isha_value;
    double isha_diff = night_portion(isha_angle) * night_time;
    if (isnan(times[Isha]) || time_diff(times[Sunset], times[Isha]) > isha_diff)
        times[Isha] = times[Sunset] + isha_diff;

    // Adjust Maghrib
    double maghrib_angle = method_params[calc_method].maghrib_is_minutes ? 4.0 : method_params[calc_method].maghrib_value;
    double maghrib_diff = night_portion(maghrib_angle) * night_time;
    if (isnan(times[Maghrib]) || time_diff(times[Sunset], times[Maghrib]) > maghrib_diff)
        times[Maghrib] = times[Sunset] + maghrib_diff;
}
/* adjust times in a prayer time array */
void adjust_times(double times[])
{
    for (int i = 0; i < TimesCount; ++i)
        times[i] += timezone - longitude / 15.0;
    times[Dhuhr] += dhuhr_minutes / 60.0;        // Dhuhr
    if (method_params[calc_method].maghrib_is_minutes)        // Maghrib
        times[Maghrib] = times[Sunset] + method_params[calc_method].maghrib_value / 60.0;
    if (method_params[calc_method].isha_is_minutes)        // Isha
        times[Isha] = times[Maghrib] + method_params[calc_method].isha_value / 60.0;

    if (adjust_high_lats != None)
        adjust_high_lat_times(times);
}
/* convert hours to day portions  */
void day_portion(double times[])
{
    for (int i = 0; i < TimesCount; ++i)
        times[i] /= 24.0;
}

/* return prayer times for a given date */
void get_prayer_times(int year, int month, int day, double _latitude, double _longitude, double _timezone, double times[])
{
    latitude = _latitude;
    longitude = _longitude;
    timezone = _timezone;
    julian_date = get_julian_date(year, month, day) - longitude / (double)(15 * 24);
    compute_day_times(times);
}

void set_calc_method(CalculationMethod method_id)
{
    calc_method = method_id;
}

/* set the angle for calculating Fajr */
void set_fajr_angle(double angle)
{
    method_params[Custom].fajr_angle = angle;
    calc_method = Custom;
}

/* set the angle for calculating Maghrib */
void set_maghrib_angle(double angle)
{
    method_params[Custom].maghrib_is_minutes = false;
    method_params[Custom].maghrib_value = angle;
    calc_method = Custom;
}

/* set the angle for calculating Isha */
void set_isha_angle(double angle)
{
    method_params[Custom].isha_is_minutes = false;
    method_params[Custom].isha_value = angle;
    calc_method = Custom;
}

/* set the minutes after mid-day for calculating Dhuhr */
void set_dhuhr_minutes(double minutes)
{
    dhuhr_minutes = minutes;
}

/* set the minutes after Sunset for calculating Maghrib */
void set_maghrib_minutes(double minutes)
{
    method_params[Custom].maghrib_is_minutes = true;
    method_params[Custom].maghrib_value = minutes;
    calc_method = Custom;
}

/* set the minutes after Maghrib for calculating Isha */
void set_isha_minutes(double minutes)
{
    method_params[Custom].isha_is_minutes = true;
    method_params[Custom].isha_value = minutes;
    calc_method = Custom;
}



std::string decimalToTime(double decimalHours) {
    // Extracting the integer part (hours) and the decimal part (minutes)
    int hours = static_cast<int>(decimalHours);
    int minutes = static_cast<int>((decimalHours - hours) * 60 + 0.5); // Adding 0.5 for rounding

    // Formatting hours and minutes
    std::string hoursStr = (hours < 10) ? "0" + std::to_string(hours) : std::to_string(hours);
    std::string minutesStr = (minutes < 10) ? "0" + std::to_string(minutes) : std::to_string(minutes);

    // Concatenating hours and minutes with a colon
    return hoursStr + ":" + minutesStr;
}

int main() {
    // Set parameters for calculation
    set_calc_method(Jafari); // Calculation method
    set_asr_method(Shafii); // Juristic method for Asr
    set_high_lats_adjust_method(AngleBased); // Adjusting method for higher latitudes

    // Set the coordinates and timezone
    double latitude = 48.8575; // Paris
    double longitude = 2.3514; // Paris
    double timezone = 1;

    // Set the date
    int year = 2024;
    int month = 5;
    int day = 2;

    // Array to hold prayer times
    double times[TimesCount];

    // Get prayer times
    get_prayer_times(year, month, day, latitude, longitude, timezone, times);

    // Display prayer times

    std::cout << "Fajr: " << decimalToTime(times[Fajr]) << " /// correct time: 04:41 or 05:00" << std::endl;
    std::cout << "Dhuhur: " << decimalToTime(times[Dhuhr]) << " /// correct time: 13:49 or 13:52" << std::endl;
    std::cout << "Asr: " << decimalToTime(times[Asr]) << " /// correct time: 17:48 or 17:49" << std::endl;
    std::cout << "Maghrib: " << decimalToTime(times[Maghrib]) << " /// correct time: 21:10 or 21:12" << std::endl;
    std::cout << "Isha: " << decimalToTime(times[Isha]) << " /// correct time: 22:40 or 22:34" << std::endl;

    std::cout << "\nSunrise: " << decimalToTime(times[Sunrise]) << " ///" << std::endl;
    std::cout << "Sunset: " << decimalToTime(times[Sunset]) << " ///" << std::endl;

    
    return 0;

}
