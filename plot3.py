import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

eta_values = np.arange(0.1, 1.0, 0.1)

colors = plt.cm.viridis(np.linspace(0, 1, len(eta_values)))

plt.figure(figsize=(10, 6))

for i, eta in enumerate(eta_values):
    filename = f"order_parameter_{eta:.1f}.txt"
    
    try:
       
        data = pd.read_csv(filename)
        time = data.iloc[:, 0]
        velocity = data.iloc[:, 1]
        plt.plot(time, velocity, color=colors[i], label=f'η = {eta:.1f}')
        
    except FileNotFoundError:
        print(f"Warning: File {filename} not found")

plt.xlabel('Time')
plt.ylabel('Average Velocity')
plt.title('Average Velocity vs Time for Different η Values')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()

plt.savefig('velocity_plots.png', bbox_inches='tight', dpi=300)

plt.show()
