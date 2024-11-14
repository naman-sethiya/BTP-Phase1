import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Initialize lists to store results
eta_values = np.arange(0.1, 1.0, 0.1)
average_order_parameters = []

# Process each file
for eta in eta_values:
    try:
        # Read the CSV file
        filename = f'order_parameter_{eta:.1f}.txt'
        data = pd.read_csv(filename)
        
        # Assuming first column is time and second is order parameter
        time = data.iloc[:, 0]
        order_parameter = data.iloc[:, 1]
        
        # Filter data for time between 5000 and 10000
        mask = (time >= 5000) & (time <= 10000)
        filtered_order_parameter = order_parameter[mask]
        
        # Calculate average order parameter for this eta
        avg_order = np.mean(filtered_order_parameter)
        average_order_parameters.append(avg_order)
        
        print(f'Eta = {eta:.1f}, Average Order Parameter = {avg_order:.4f}')
        
    except FileNotFoundError:
        print(f"Warning: File {filename} not found")
    except Exception as e:
        print(f"Error processing eta = {eta:.1f}: {str(e)}")
        average_order_parameters.append(np.nan)

# Convert to numpy array for easier handling
average_order_parameters = np.array(average_order_parameters)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot without error bars
plt.plot(eta_values, average_order_parameters, 'bo-', linewidth=2, markersize=8)

# Customize the plot
plt.xlabel('Noise (η)', fontsize=12)
plt.ylabel('Order Parameter (φ)', fontsize=12)
plt.title('Order Parameter vs Noise', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)

# Set axis limits with some padding
plt.xlim(0.05, 0.95)
y_min = np.min(average_order_parameters) * 0.9
y_max = np.max(average_order_parameters) * 1.1
plt.ylim(y_min, y_max)

# Save the plot
plt.savefig('average_order_parameter_vs_eta.png', dpi=300, bbox_inches='tight')
plt.show()

# Save the numerical data to a CSV file
results_df = pd.DataFrame({
    'eta': eta_values,
    'phi': average_order_parameters
})
results_df.to_csv('average_order_parameters.csv', index=False)

