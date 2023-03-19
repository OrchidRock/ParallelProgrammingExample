#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <gmpxx.h>
#include <openmpi/mpi.h>

using namespace std;

#define ll long int
#define DIMENSION 2     // The dimension of Newton-Space.
#define BUFFER_MAX_LEN 4096

const mpf_class G("6.673e-11"); //  The gravitational constant: 6.673e-11

enum Coordinate {
    X=0, Y=1, Z=2
};

struct Vector {
    mpf_class x;
    mpf_class y;
    
    Vector(): x(0), y(0) {}
    Vector(ll xx, ll yy) : x(xx), y(yy) {}
    Vector(string xx, string yy) : x(xx), y(yy) {}
    Vector(mpf_class xx, mpf_class yy) : x(xx), y(yy) {}
    
    Vector operator-(const Vector &b) {
        return { x - b.x, y - b.y };
    }
    Vector operator+(const Vector &b) {
        return { x + b.x, y + b.y };
    }
    Vector operator*(ll scalar) {
        return { x*scalar, y*scalar };
    }
    Vector operator*(mpf_class scalar) {
        return { x*scalar, y*scalar };
    }
    
    void operator+=(const Vector &b) {
        x = x + b.x;
        y = y + b.y;
    }
    void operator-=(const Vector &b) {
        x = x - b.x;
        y = y - b.y;
    }
    void operator*=(const Vector &b) {
        x = x * b.x;
        y = y * b.y;
    }
    void operator/=(const Vector &b) {
        x = x / b.x;
        y = y / b.y;
    }
    
    mpf_class value() {
        return sqrt(x*x + y*y);
    }
};

// The State of Particle.
struct ParticleState {
    Vector pos;    // the position
    Vector vel;    // the velocity
    Vector force;   // the total gravitation force
};


void get_input(int my_rank, char *buffer, int *buffer_len_ptr) {

    if(my_rank == 0) {
        *buffer_len_ptr = fread(buffer, 1, BUFFER_MAX_LEN, stdin);
    }
    MPI_Bcast(buffer, BUFFER_MAX_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(buffer_len_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);


    // Get Input Data
    int comm_sz, my_rank;
    int n = 0;  // The number of particles
    int delta;  // The timestep
    int T;      // The number of timesteps
    
    char buffer[BUFFER_MAX_LEN];
    int buffer_len_ptr = 0;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

    get_input(my_rank, buffer, &buffer_len_ptr);
    
    stringstream str_cin(buffer); 
    str_cin >> n >> delta >> T;

    vector<mpf_class> masses(n, 0); // the mass of particle.
    vector<ParticleState> pstates(n);
    string para;
    for(int i=0; i<n; i++){
        str_cin >> para;
        masses[i] = para;
    }

    string para_x, para_y;
    // The initial state
    for(int i = 0; i < n; i++) {
        str_cin >> para_x >> para_y;

        pstates[i].pos = {para_x, para_y};
        
        str_cin >> para_x >> para_y;
        pstates[i].vel = {para_x, para_y};

        pstates[i].force  = {0, 0};
    }
    

    mpf_class dist, dist_cubed;
    Vector diff(0,0);
    Vector force_jk(0,0);
    
    int local_n = n / comm_sz;
    int local_a = local_n*my_rank;
    int local_b = (my_rank != (n-1)) ? (local_a + local_n) : n;
    

    for(int i = 0; i < T; i++) {

        // Compute total force to each particle j.
        for(int j = local_a; j < local_b; j++) {
            pstates[j].force = {0, 0};
        }
        for(int j = local_a; j < local_b; j++) {
            for(int k = j + 1; k < n; k++){
                
                diff = pstates[k].pos - pstates[j].pos; 
                dist = diff.value();
                dist_cubed = dist*dist*dist;
                force_jk = diff * ((G * masses[j] * masses[k]) / dist_cubed);
                
                // Use Newton's  third law of motion.
                pstates[j].force += force_jk;

                pstates[k].force -= force_jk; 
            }
        }
        
        //MPI_Allgather()
         
        /* 
        cout << "Force: " << endl;
        for(int j = 0; j < n; j++) {
            cout << pstates[j].force.x << " " << pstates[j].force.y << "\n";
        }
        cout << endl;    
        */

        // Compute position and velocity of each particle j.
        // Euler's Method
        for(int j = local_a; j < local_b; j++) {
            pstates[j].pos += pstates[j].vel * delta;
            pstates[j].vel += pstates[j].force * (delta / masses[j]);
        }
        
    }
    
    // outout
    if(my_rank == 0) {
        cout << "Final State: " << endl; 
        for(int j = 0; j < n; j++) {
            cout << pstates[j].pos.x << " " << pstates[j].pos.y << " ";
            cout << pstates[j].vel.x << " " << pstates[j].vel.y << "\n";
        }
        cout << endl;    
    }
    MPI_Finalize();
    return 0;
}
