module Time

using Dates: Date, DateTime, format, isleapyear, now

# global constants
const MDAYS_LEAP  = [0,31,60,91,121,152,182,213,244,274,305,335,366];
const MDAYS       = [0,31,59,90,120,151,181,212,243,273,304,334,365];
const NDAYS_LEAP  = [31,29,31,30,31,30,31,31,30,31,30,31];
const NDAYS       = [31,28,31,30,31,30,31,31,30,31,30,31];
const TIME_FORMAT = ["YYYYMMDD", "YYYYMMDDhh", "YYYYMMDDhhmm", "YYYYMMDDhhmmss"];
const TIME_OUTPUT = ["DATE", "DATETIME", "DOY", "FDOY"];


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Aug-24: move function outside of the folder
#     2022-Aug-24: add method to convert float doy
#
#######################################################################################################################################################################################################
"""

    parse_timestamp(timestamp::Union{Int,String}; in_format::String = "YYYYMMDD", out_format::String = "DOY")
    parse_timestamp(year::Int, doy::Int; out_format::String = "DOY")
    parse_timestamp(year::Int, doy::AbstractFloat; out_format::String = "DOY")

Convert timestamp, given
- `timestamp` Time stamp
- `in_format` Format of `timestamp`, default is `YYYYMMDD`
- `out_format` Output format, default is `DOY`
- `year` Year (in this case, the function will convert year and day to timestamp first)
- `doy` Day of year (typically 1-365, 1-366 for leap years)

The input format (string or integer) supports `YYYYMMDD`, `YYYYMMDDhh`, `YYYYMMDDhhmm`, and `YYYYMMDDhhmmss`, where the labels are
- `YYYY` Year number
- `MM` Month number
- `DD` Day number
- `hh` Hour number
- `mm` Minute number
- `ss` second number

The supported outputs are
- `DATE` A `Dates.Date` type variable
- `DATETIME` A `Dates.DateTime` type variable
- `DOY` A day of year integer
- `FDOY` A day of year float

---
Examples
```julia
time = parse_timestamp(20200130; in_format="YYYYMMDD", out_format="FDOY");
time = parse_timestamp("20200130"; in_format="YYYYMMDD", out_format="FDOY");
time = parse_timestamp(2020, 100);
time = parse_timestamp(2020, 100.23435436);
```

"""
function parse_timestamp end

parse_timestamp(timestamp::Union{Int,String}; in_format::String = "YYYYMMDD", out_format::String = "DOY") = (
    _time_stamp = string(timestamp);
    @assert in_format in TIME_FORMAT "in_format not found in PkgUtility.TIME_FORMAT";
    @assert out_format in TIME_OUTPUT "out_format not found in PkgUtility.TIME_OUTPUT";

    # get year, month, and day
    _year   = parse(Int, _time_stamp[1:4]);
    _month  = parse(Int, _time_stamp[5:6]);
    _day    = parse(Int, _time_stamp[7:8]);
    _hour   = 0;
    _minute = 0;
    _second = 0;
    @assert 0 <= _year  <= 9999;
    @assert 1 <= _month <= 12;
    @assert 1 <= _day   <= month_days(_year,_month);

    # get hour, minute, and second
    if in_format in ["YYYYMMDDhh", "YYYYMMDDhhmm", "YYYYMMDDhhmmss"]
        _hour   = parse(Int, _time_stamp[9:10]);
        @assert 0 <= _hour <= 24;
    end;

    if in_format in ["YYYYMMDDhhmm", "YYYYMMDDhhmmss"]
        _minute = parse(Int, _time_stamp[11:12]);
        @assert 0 <= _minute <= 59;
    end;

    if in_format in ["YYYYMMDDhhmmss"]
        _second = parse(Int, _time_stamp[13:14]);
        @assert 0 <= _second <= 59;
    end;

    # determine the format to return (float ddoy by default)
    if out_format=="DATE"
        return Date(_year, _month, _day)
    end;

    if out_format=="DATETIME"
        return DateTime(_year, _month, _day, _hour, _minute, _second)
    end;

    _doy = (isleapyear(_year) ? MDAYS_LEAP[_month] : MDAYS[_month]) + _day;

    if out_format=="DOY"
        return _doy
    end;

    return _doy + (_hour + _minute/60 + _second/3600) / 24
);

parse_timestamp(year::Int, doy::Int; out_format::String = "DOY") = (
    @assert 0 <= year <= 9999;
    @assert 1 <= doy  <= (isleapyear(year) ? 366 : 365);

    # convert the year and doy to timestamp
    _month = month_ind(year, doy);
    _day   = doy - (isleapyear(year) ? MDAYS_LEAP[_month] : MDAYS[_month]);
    _stamp = lpad(year, 4, "0") * lpad(_month, 2, "0") * lpad(_day, 2, "0");

    return parse_timestamp(_stamp; in_format = "YYYYMMDD", out_format = out_format)
);

parse_timestamp(year::Int, doy::AbstractFloat; out_format::String = "DOY") = (
    @assert 0 <= year <= 9999;
    @assert 1 <= doy  <  (isleapyear(year) ? 367 : 366);

    # convert the year and doy to timestamp
    _doy   = Int(floor(doy));
    _month = month_ind(year, _doy);
    _day   = _doy - (isleapyear(year) ? MDAYS_LEAP[_month] : MDAYS[_month]);
    _fhour = (doy % 1) * 24;
    _hour  = Int(floor(_fhour));
    _fminu = (_fhour % 1) * 60;
    _minu  = Int(floor(_fminu));
    _fsec  = (_fminu % 1) * 60;
    _sec   = Int(floor(_fsec));
    _stamp = lpad(year, 4, "0") * lpad(_month, 2, "0") * lpad(_day, 2, "0") * lpad(_hour, 2, "0") * lpad(_minu, 2, "0") * lpad(_sec, 2, "0");

    return parse_timestamp(_stamp; in_format = "YYYYMMDDhhmmss", out_format = out_format)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Aug-24: move function outside of the folder
#
#######################################################################################################################################################################################################
"""

    month_days(year::Int, month::Int)

Return the number of days per month, given
- `year` Year
- `month` Month

"""
function month_days(year::Int, month::Int)
    @assert 1 <= month <= 12;

    return isleapyear(year) ? NDAYS_LEAP[month] : NDAYS[month]
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Aug-24: move function outside of the folder
#     2022-Aug-24: add method to convert float doy
#
#######################################################################################################################################################################################################
"""

    month_ind(year::Int, doy::Int)
    month_ind(year::Int, doy::AbstractFloat)

Return the month index, given
- `year` Year
- `doy` Day of year (typically 1-365, 1-366 for leap years)

"""
function month_ind end

month_ind(year::Int, doy::Int) = (
    # if is leap year
    if isleapyear(year)
        @assert 1 <= doy <= 366;
        _month = 1;
        for i in 1:12
            if MDAYS_LEAP[i] < doy <= MDAYS_LEAP[i+1]
                _month = i;
                break;
            end;
        end;

        return _month
    end;

    # if not leap year
    @assert 1 <= doy <= 365;
    _month = 1;
    for i in 1:12
        if MDAYS[i] < doy <= MDAYS[i+1]
            _month = i;
            break;
        end;
    end;

    return _month
);

month_ind(year::Int, doy::AbstractFloat) = month_ind(year, Int(floor(doy)));


end # module
