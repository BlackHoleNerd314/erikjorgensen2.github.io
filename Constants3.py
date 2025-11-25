import numpy as np
#
tU = 4.35456 * 10**17
c = 2.99792458 * 10**8
G = 6.67408 / 10**11
h = 6.62607015 / 10**34
kB = 1.380649 / 10**23
Qe = 1.602176634 / 10**19
Ke = (c**2) / 10**7
Me = (5.485799090441 / 10**4) / (6.02214076 * 10**26)
#
tU0 = [0,0,1,0,0]
c0 = [1,0,-1,0,0]
G0 = [3,-1,-2,0,0]
h0 = [2,1,-1,0,0]
kB0 = [2,1,-2,-1,0]
Qe0 = [0,0,0,0,1]
Ke0 = [3,1,-2,0,-2]
Me0 = [0,1,0,0,0]
#
c1 = np.log10(c)
G1 = np.log10(G)
h1 = np.log10(h)
Gh1 = G1 + h1
tU1 = np.log10(tU)

def is_leap_year(year):
    # Gregorian calendar leap year rules
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

def days_in_month(year, month):
    if month == 2:
        return 29 if is_leap_year(year) else 28
    return [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31][month - 1]

def weekday_name(index):
    # 0 = Monday, 6 = Sunday (ISO 8601), but we’ll map Saturday as base
    return ["Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday"][index % 7]

def convert_seconds_to_gregorian(total_seconds):
    SECONDS_PER_DAY = 86400

    # Split total seconds
    total_days = total_seconds // SECONDS_PER_DAY
    seconds_in_day = total_seconds % SECONDS_PER_DAY

    # Compute weekday assuming 0000-01-01 is Saturday
    weekday = weekday_name(total_days % 7)

    # Date calculation
    year = 0
    days_remaining = total_days
    while True:
        days_in_year = 366 if is_leap_year(year) else 365
        if days_remaining >= days_in_year:
            days_remaining -= days_in_year
            year += 1
        else:
            break

    month = 1
    while True:
        dim = days_in_month(year, month)
        if days_remaining >= dim:
            days_remaining -= dim
            month += 1
        else:
            break

    day = days_remaining + 1

    hour = seconds_in_day // 3600
    minute = (seconds_in_day % 3600) // 60
    second = seconds_in_day % 60

    return f"{weekday}, {year:04d}-{month:02d}-{day:02d} {hour:02d}:{minute:02d}:{second:02d}"

import requests

def reverse_geocode(lat, lon):
    url = "https://nominatim.openstreetmap.org/reverse"
    params = {
        "lat": lat,
        "lon": lon,
        "format": "jsonv2",
        "addressdetails": 1
    }
    headers = {
        "User-Agent": "YourApp/1.0 (email@example.com)"
    }
    response = requests.get(url, params=params, headers=headers)
    if response.ok:
        data = response.json()
        return {
            "display_name": data.get("display_name", "Unknown"),
            "address": data.get("address", {}),
            "lat": data.get("lat"),
            "lon": data.get("lon")
        }
    else:
        return {"error": f"Failed with code {response.status_code}"}


# Example 1: Land
print(reverse_geocode(38.8568, -77.2345))
print(convert_seconds_to_gregorian(63115608251))

N = 0
N0 = 201
time0 = np.sqrt(h*G/c**5)
space0 = np.sqrt(h*G/c**3)

for time in range(0,N0):
    for space in range(0,N0):
        for matter in range(-N0,N0):
            c2 = time - space
            G2 = - 2*time + 3*space - matter
            h2 = - time + 2*space + matter
            if (c2>=0)&(G2>=0)&(h2>=0):
                print("Time:" + str(time))
                print("Space:" + str(space))
                print("Matter:" + str(matter))

                N = N + 1

print(N)