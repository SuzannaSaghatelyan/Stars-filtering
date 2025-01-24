import csv
from datetime import datetime
from math import sqrt, radians, cos, sin, atan2


def parse_tsv(file_path):
    """Parse TSV file and return its data as a list of dictionaries"""
    with open(file_path, 'r') as tsv_file:
        lines = tsv_file.readlines()

        data_lines = [line for line in lines if not line.startswith('#')]

        if data_lines:
            reader = csv.DictReader(data_lines, delimiter='\t')

            return [row for row in reader]
        else:
            raise ValueError("No valid data found in the file.")


def haversine_distance(ra1, dec1, ra2, dec2):
    """Calculate angular distance between two points on a sphere using RA and DEC"""
    ra1, dec1, ra2, dec2 = map(radians, [ra1, dec1, ra2, dec2])
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    a = sin(delta_dec / 2) ** 2 + cos(dec1) * cos(dec2) * sin(delta_ra / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return c

def quick_sort(arr, key_func):
    """Implement quick sort algorithm"""
    if len(arr) <= 1:
        return arr
    pivot = arr[len(arr) // 2]
    pivot_key = key_func(pivot)
    left = [x for x in arr if key_func(x) < pivot_key]
    middle = [x for x in arr if key_func(x) == pivot_key]
    right = [x for x in arr if key_func(x) > pivot_key]
    return quick_sort(left, key_func) + middle + quick_sort(right, key_func)

def filter_stars(data, ra, dec, fov_h, fov_v):
    """Filter stars within the fov defined by ra, dec, fov_h, and fov_v"""
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

def find_brightest_stars(stars, ra, dec, n):
    """Find the n brightest stars sorted by distance from the given point"""
    for star in stars:
        try:
            star['magnitude'] = float(star['phot_g_mean_mag'])
            star['distance'] = haversine_distance(ra, dec, float(star['ra_ep2000']), float(star['dec_ep2000']))
        except ValueError:
            star['magnitude'] = float('inf')
            star['distance'] = float('inf')

    stars = quick_sort(stars, key_func=lambda x: (x['distance'], x['magnitude']))
    return stars[:n]

def save_to_csv(stars, output_file):
    """Save the stars to CSV file"""
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

def main():
    ra = float(input("Enter RA "))
    dec = float(input("Enter DEC "))
    fov_h = float(input("Enter FOV horizontal "))
    fov_v = float(input("Enter FOV vertical "))
    n = int(input("Enter number of stars "))

    data = parse_tsv("data/stars_small.tsv")

    filtered_stars = filter_stars(data, ra, dec, fov_h, fov_v)

    brightest_stars = find_brightest_stars(filtered_stars, ra, dec, n)

    output_file = datetime.now().strftime("%Y%m%d_%H%M%S.csv")
    save_to_csv(brightest_stars, output_file)

if __name__ == "__main__":
    main()
