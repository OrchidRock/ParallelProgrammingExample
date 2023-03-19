#include <iostream>
#include <vector>
#include <cmath>
#include <gmpxx.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define ll long int
#define DIMENSION 2     // The dimension of Newton-Space.

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

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    // Get Input Data

    int n = 0;      // The number of particles
    int delta_t;    // The timestep
    int T;          // The number of timesteps
    string para;    // the argument

    cin >> n >> delta_t >> T;
    
    vector<mpf_class> masses(n, 0); // the mass of particle.
    vector<ParticleState> pstates(n);

#ifdef _OPENMP
    vector<omp_lock_t> locks(n);
    for(int i=0; i<n; i++) {
        omp_init_lock(&locks[i]);
    }
#endif

    for(int i = 0 ; i < n; i++) {
        cin >> para;
        masses[i] = para;
    }
    

    string para_x, para_y;
    // The initial state
    for(int i = 0; i < n; i++) {
        cin >> para_x >> para_y;

        pstates[i].pos = {para_x, para_y};
        
        cin >> para_x >> para_y;
        pstates[i].vel = {para_x, para_y};

        pstates[i].force  = {0, 0};
    }
    

    mpf_class dist, dist_cubed;
    Vector diff(0,0);
    Vector force_jk(0,0);

    #pragma omp parallel num_threads(n)
    for(int i = 0; i < T; i++) {

        // Compute total force to each particle j.
        #pragma omp for
        for(int j = 0; j < n; j++) {
            pstates[j].force = {0, 0};
        }

        #pragma omp for private(diff, dist, dist_cubed, force_jk) schedule(dynamic)
        for(int j = 0; j < n; j++) {
            for(int k = j + 1; k < n; k++){
                
                //cout << "j: "<<j  << " k: " << k << endl; 
                diff = pstates[k].pos - pstates[j].pos; 
                dist = diff.value();
                dist_cubed = dist*dist*dist;
                force_jk = diff * ((G * masses[j] * masses[k]) / dist_cubed);
                
                // Use Newton's  third law of motion.

#ifdef _OPENMP
                omp_set_lock(&locks[j]);
                pstates[j].force += force_jk;
                omp_unset_lock(&locks[j]);

                omp_set_lock(&locks[k]);
                pstates[k].force -= force_jk;
                omp_unset_lock(&locks[k]);
#else
                #pragma omp critical
                {
                    pstates[j].force += force_jk;
                    pstates[k].force -= force_jk; 
                }
#endif

            }
        }
        
       /* 
        cout << "Force: " << endl;
        for(int j = 0; j < n; j++) {
            cout << pstates[j].force.x << " " << pstates[j].force.y << "\n";
        }
        cout << endl;    
        */

        // Compute position and velocity of each particle j.
        // Euler's Method
        #pragma omp for
        for(int j = 0; j < n; j++) {
            pstates[j].pos += pstates[j].vel * delta_t;
            pstates[j].vel += pstates[j].force * (delta_t / masses[j]);
        }
        
    }
    
    // outout
    cout << "Final State: " << endl; 
    for(int j = 0; j < n; j++) {
        cout << pstates[j].pos.x << " " << pstates[j].pos.y << " ";
        cout << pstates[j].vel.x << " " << pstates[j].vel.y << "\n";
    }
    cout << endl;

#ifdef _OPENMP
    for(int i=0; i<n; i++) {
        omp_destroy_lock(&locks[i]);
    }
#endif
    return 0;
}
