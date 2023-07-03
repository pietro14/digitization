import csv

# the electric field is in kV/cm
# the drift velocity is in cm/microseconds
def compute_drift_velocity(electric_field):
    with open('drift_velocity_datasets/drift_velocity_HeCF4_60_40.csv', 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  # Skip the header
        data = list(reader)
        
        # Find the closest lower and upper electric field values
        lower_field, upper_field = None, None
        for row in data:
            field = float(row[0])
            if field < electric_field:
                lower_field = field
            else:
                upper_field = field
                break
        
        # If there is no lower or upper field, return None
        if lower_field is None or upper_field is None:
            return None
        
        # Find the corresponding drift velocities
        lower_velocity, upper_velocity = None, None
        for row in data:
            field = float(row[0])
            velocity = float(row[1])
            if field == lower_field:
                lower_velocity = velocity
            elif field == upper_field:
                upper_velocity = velocity
        
        # If there is no lower or upper velocity, return None
        if lower_velocity is None or upper_velocity is None:
            return None
        
        # Perform linear interpolation
        drift_velocity = lower_velocity + (electric_field - lower_field) * \
                         (upper_velocity - lower_velocity) / (upper_field - lower_field)
        
        return drift_velocity

