import csv

def compute_value_from_csv(file_path, x, x_column, value_column):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        data = list(reader)

        # Find the closest lower and upper x values
        lower_x, upper_x = None, None
        for row in data:
            curr_x = float(row[x_column])
            if curr_x < x:
                lower_x = curr_x
            else:
                upper_x = curr_x
                break

        # If there is no lower or upper x, return None
        if lower_x is None or upper_x is None:
            return None

        # Find the corresponding value
        lower_value, upper_value = None, None
        for row in data:
            curr_x = float(row[x_column])
            value = float(row[value_column])
            if curr_x == lower_x:
                lower_value = value
            elif curr_x == upper_x:
                upper_value = value

        # If there is no lower or upper value, return None
        if lower_value is None or upper_value is None:
            return None

        # Perform linear interpolation
        interpolated_value = lower_value + (x - lower_x) * \
                             (upper_value - lower_value) / (upper_x - lower_x)

        return interpolated_value

