import csv
from datetime import datetime
from math import sqrt, radians, cos, sin, atan2

def is_valid_ra_dec(ra, dec):
    if not (0 <= ra <= 360):
        return False
    if not (-90 <= dec <= 90):
        return False
    return True

def read_tsv_file(file_path: str) -> list[dict[str, str]]:
    with open(file_path, 'r') as tsv_file:
        next(tsv_file)
        rows = csv.DictReader(tsv_file, delimiter='\t')

        stars_data = []
        for row in rows:
            try:
                ra = float(row['ra_ep2000'].strip())
                dec = float(row['dec_ep2000'].strip())

                if not is_valid_ra_dec(ra, dec):
                    continue

                star = {
                    'source_id': row['source_id'].strip(),
                    'ra_ep2000': ra,
                    'dec_ep2000': dec,
                    'brightness': float(row['phot_g_mean_mag'].strip()),
                }
                stars_data.append(star)
            except (ValueError, KeyError):
                continue

        return stars_data


def calculate_haversine_distance(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(radians, [ra1, dec1, ra2, dec2])
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    a = sin(delta_dec / 2) ** 2 + cos(dec1) * cos(dec2) * sin(delta_ra / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return c

def quick_sort(arr, key_func):
    if len(arr) <= 1:
        return arr
    pivot = arr[len(arr) // 2]
    pivot_key = key_func(pivot)
    left = [x for x in arr if key_func(x) < pivot_key]
    middle = [x for x in arr if key_func(x) == pivot_key]
    right = [x for x in arr if key_func(x) > pivot_key]
    return quick_sort(left, key_func) + middle + quick_sort(right, key_func)


def filter_stars_in_fov(data, ra, dec, fov_h, fov_v):
    filtered = []
    for star in data:
        try:
            star_ra = float(star['ra_ep2000'])
            star_dec = float(star['dec_ep2000'])
            if abs(star_ra - ra) <= fov_h / 2 and abs(star_dec - dec) <= fov_v / 2:
                filtered.append(star)
        except ValueError:
            continue
    return filtered

def filter_stars_by_distance(stars, ra, dec):
    for star in stars:
        try:
            star['distance'] = calculate_haversine_distance(ra, dec, float(star['ra_ep2000']), float(star['dec_ep2000']))
        except ValueError:
            star['distance'] = float('inf')

    stars = quick_sort(stars, key_func=lambda x: x['distance'])
    return stars


def sort_stars_by_brightness(stars):
    for star in stars:
        try:
            star['magnitude'] = float(star['brightness'])
        except ValueError:
            star['magnitude'] = float('inf')

    stars = quick_sort(stars, key_func=lambda x: x['magnitude'])
    return stars

def save_to_csv(stars):
    output_file = datetime.now().strftime("%Y%m%d_%H%M%S.csv")

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["ID", "RA", "DEC", "Brightness(Magnitude)", "Distance"])
        for star in stars:
            writer.writerow([
                star['source_id'],
                star['ra_ep2000'],
                star['dec_ep2000'],
                star['magnitude'],
                star['distance']
            ])

def get_valid_input(prompt, min_value, max_value):
    while True:
        try:
            value = float(input(prompt))
            if min_value <= value <= max_value:
                return value
            else:
                print(f"Input out of range. Please enter a value between {min_value} and {max_value}")
        except ValueError:
            print("Invalid input.")

def main():
    ra = get_valid_input("Enter RA (0-360) ", 0, 360)
    dec = get_valid_input("Enter DEC (-90 to 90) ", -90, 90)
    fov_h = float(input("Enter FOV horizontal "))
    fov_v = float(input("Enter FOV vertical "))
    n = int(input("Enter number of stars to find: "))

    data = read_tsv_file("data/stars_small.tsv")

    filtered_stars = filter_stars_in_fov(data, ra, dec, fov_h, fov_v)

    stars_by_distance = filter_stars_by_distance(filtered_stars, ra, dec)

    brightest_stars = sort_stars_by_brightness(stars_by_distance)

    brightest_stars = brightest_stars[:n]
    save_to_csv(brightest_stars)

if __name__ == "__main__":
    main()
