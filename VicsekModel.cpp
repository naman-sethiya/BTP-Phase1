#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <array>
#include <fstream>

#define PI 3.14159265358979323846
#define N 1024          // Number of particles
#define L 32.0         // Size of the box
#define R 1.0           // Interaction radius
#define V0 1.0          // Speed of particles
#define TIMESTEPS 20000 // Number of time steps
#define DT 1.0          // Time step size
#define ETA 0.3        // Noise strength

// Cell list related constants
#define CELL_SIZE R     // Cell size equals interaction radius
#define NUM_CELLS (int)(L/CELL_SIZE) // Number of cells in each dimension

// RNG constants (unchanged)
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

using namespace std;

struct Particle {
    float x, y;    // Position
    float theta;   // Direction
    float vx, vy;  // Velocity components
};

// Cell list data structure
vector<int> head;  
vector<int> list;  

float ran2(long *idum) {
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {  // Initialize
        if (-(*idum) < 1) *idum = 1;  // Prevent idum = 0
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; j--) {  // Load the shuffle table (after 8 warm-ups)
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;  // Start here when not initializing
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

// Get cell index from position
int get_cell_index(float x, float y) {
    int cx = int(floor(x / CELL_SIZE)) % NUM_CELLS;
    int cy = int(floor(y / CELL_SIZE)) % NUM_CELLS;
    
    if(cx < 0) cx += NUM_CELLS;
    if(cy < 0) cy += NUM_CELLS;
    
    return cy * NUM_CELLS + cx;
}

// Build cell lists
void build_cell_lists(const vector<Particle>& particles) {
    fill(head.begin(), head.end(), -1);
    
    for(int i = 0; i < N; i++) {
        int cell = get_cell_index(particles[i].x, particles[i].y);
        list[i] = head[cell];
        head[cell] = i;
    }
}

// Calculate order parameter 
float calculate_order_parameter(const vector<Particle>& particles) {
    float vx_sum = 0, vy_sum = 0;
    for(const auto& p : particles) {
        vx_sum += V0 * cos(p.theta);
        vy_sum += V0 * sin(p.theta);
    }
    float magnitude = sqrt(vx_sum*vx_sum + vy_sum*vy_sum);
    return magnitude / (N * V0); // Normalize by N*V0
}

// Initialize particles with more randomness
void initialize_particles(vector<Particle>& particles, long *idum) {
    for(int i = 0; i < N; i++) {
        particles[i].x = L * ran2(idum);
        particles[i].y = L * ran2(idum);
        // Initialize with completely random directions
        particles[i].theta = 2 * PI * ran2(idum);
        particles[i].vx = V0 * cos(particles[i].theta);
        particles[i].vy = V0 * sin(particles[i].theta);
    }
}

// Periodic boundary distance calculation
float periodic_distance(float dx, float dy) {
    if(dx > L/2) dx -= L;
    else if(dx < -L/2) dx += L;
    if(dy > L/2) dy -= L;
    else if(dy < -L/2) dy += L;
    return sqrt(dx*dx + dy*dy);
}

// Calculate average direction using cell lists 
float average_direction(const Particle& p, const vector<Particle>& particles, long* idum) {
    float sum_sin = 0.0, sum_cos = 0.0;
    int neighbor_count = 0;
    
    // Include self in the average
    sum_sin += sin(p.theta);
    sum_cos += cos(p.theta);
    neighbor_count++;
    
    // Get current cell
    int cell = get_cell_index(p.x, p.y);
    int cell_x = cell % NUM_CELLS;
    int cell_y = cell / NUM_CELLS;
    
    // Check neighboring cells (including current cell)
    for(int dy = -1; dy <= 1; dy++) {
        for(int dx = -1; dx <= 1; dx++) {
            int nx = (cell_x + dx + NUM_CELLS) % NUM_CELLS;
            int ny = (cell_y + dy + NUM_CELLS) % NUM_CELLS;
            int ncell = ny * NUM_CELLS + nx;
            
            // Traverse linked list in this cell
            for(int i = head[ncell]; i != -1; i = list[i]) {
                float dx = particles[i].x - p.x;
                float dy = particles[i].y - p.y;
                
                // Use improved periodic boundary conditions
                float dist = periodic_distance(dx, dy);
                
                if(dist < R && dist > 0) { // Exclude self from neighbor check
                    sum_sin += sin(particles[i].theta);
                    sum_cos += cos(particles[i].theta);
                    neighbor_count++;
                }
            }
        }
    }
    
    if(neighbor_count > 0) {
        // Calculate average angle correctly
        float avg_theta = atan2(sum_sin, sum_cos);
        // Add noise AFTER averaging
        float noise = ETA * (2.0 * ran2(idum) - 1.0) * PI;
        return fmod(avg_theta + noise + 2*PI, 2*PI);
    }
    
    // If no neighbors, add noise to current direction
    return fmod(p.theta + ETA * (2.0 * ran2(idum) - 1.0) * PI + 2*PI, 2*PI);
}

// Update particles with corrected timestep handling
void update_particles(vector<Particle>& particles, long *idum) {
    build_cell_lists(particles);
    
    vector<Particle> new_particles = particles;
    
    // First update all directions
    for(int i = 0; i < N; i++) {
        new_particles[i].theta = average_direction(particles[i], particles, idum);
    }
    
    // Then update all positions
    for(int i = 0; i < N; i++) {
        // Update velocity components
        new_particles[i].vx = V0 * cos(new_particles[i].theta);
        new_particles[i].vy = V0 * sin(new_particles[i].theta);
        
        // Update position
        new_particles[i].x += new_particles[i].vx * DT;
        new_particles[i].y += new_particles[i].vy * DT;
        
        // Periodic boundary conditions
        new_particles[i].x = fmod(new_particles[i].x + L, L);
        new_particles[i].y = fmod(new_particles[i].y + L, L);
    }
    
    particles = new_particles;
}

int main() {
    vector<Particle> particles(N);
    head.resize(NUM_CELLS * NUM_CELLS, -1);
    list.resize(N);
    
    long seed = -999999999;
    
    // Initialize particles
    initialize_particles(particles, &seed);
    
    // Open files for output
    ofstream order_file("order_parameter_0.3_20000.txt");
    
    // Main simulation loop
    for(int t = 0; t <= TIMESTEPS; t++) {
        update_particles(particles, &seed);
        
        if(t < 1000 && t % 10 == 0) {
            float order = calculate_order_parameter(particles);
            order_file << t << "," << order << endl;
        }
        else if(t % 100 == 0) {
            float order = calculate_order_parameter(particles);
            order_file << t << "," << order << endl;
        }
    }
    
    // Save final positions
    ofstream pos_file("final_positions_0.3_20000.txt");
    for(const auto& p : particles) {
        pos_file << p.x << "," << p.y << "," << p.theta << endl;
    }
    
    order_file.close();
    pos_file.close();
    
    return 0;
}