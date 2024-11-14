import matplotlib.pyplot as plt
import numpy as np
import csv

def read_csv(filename):
    x = []
    y = []
    theta = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:  
                x.append(float(row[0]))
                y.append(float(row[1]))
                theta.append(float(row[2]))
    return np.array(x), np.array(y), np.array(theta)

def plot_arrows(x, y, theta, box_size=32):
    plt.figure(figsize=(10, 10))
    
    arrow_scale = box_size / 100  
    head_width = box_size / 200
    head_length = box_size / 150
    
    for i in range(len(x)):
        dx = np.cos(theta[i]) * arrow_scale
        dy = np.sin(theta[i]) * arrow_scale
        plt.arrow(x[i], y[i], dx, dy, 
                 head_width=head_width, 
                 head_length=head_length, 
                 fc='blue', 
                 ec='blue',
                 alpha=0.6)  
    
    plt.xlim(0, box_size)
    plt.ylim(0, box_size)
    
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('$\eta$ = 0.9, t = 10000')
    
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    plt.axhline(y=box_size, color='k', linestyle='-', alpha=0.3)
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    plt.axvline(x=box_size, color='k', linestyle='-', alpha=0.3)
    
    plt.show()

x, y, theta = read_csv(r'C:\Users\user\OneDrive\Desktop\Sem-7\BTP\final_positions_0.9.txt')
plot_arrows(x, y, theta)
