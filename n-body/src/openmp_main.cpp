#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#define ll long long
#define DIMENSION 2     // The dimension of Newton-Space.

const ll G = 6673 /*e-11*/; //  The gravitational constant: 6.673e-11

enum Coordinate {
    X=0, Y=1, Z=2
};

struct Vector {
    ll x;
    ll y;
    Vector(): x(0), y(0) {}
    Vector(ll x, ll y): x(x), y(y) {}

    Vector operator-(const Vector &b) {
        return { x - b.x, y - b.y };
    }
    Vector operator+(const Vector &b) {
        return { x + b.x, y + b.y };
    }
    Vector operator*(ll scalar) {
        return { scalar*x, scalar*y };
    }

    double value() const{
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
    
    cin >> n >> delta_t >> T;
    
    vector<ll> masses(n, 0); // the mass of particle.
    vector<ParticleState> pstates(n);

    for(int i = 0 ; i < n; i++) {
        cin >> masses[i]; 
    }
    
    int x, y;
    // The initial state
    for(int i = 0; i < n; i++) {
        cin >> x >> y;
        pstates[i].pos = {x, y};  
        
        cin >> x >> y;
        pstates[i].vel = {x, y};

        pstates[i].force = {0, 0};
    }
    

    double dist, dist_cubed;
    Vector diff(0,0);
    Vector force_jk(0,0);

    for(int i = 0; i < T; i++) {

        // Compute total force to each particle j.
        for(int j = 0; j < n; j++) {
            pstates[j].force = {0, 0};
        }
        for(int j = 0; j < n; j++) {
            for(int k = j + 1; k < n; k++){
                
                diff = pstates[j].pos - pstates[k].pos; 
                dist = diff.value();
                dist_cubed = dist*dist*dist;
                force_jk = diff * ((G * masses[j] * masses[k]) / dist_cubed);
                
                // Use Newton's  third law of motion.
                pstates[j].force = pstates[j].force + force_jk;
                pstates[k].force = pstates[k].force - force_jk; 
            }
        }
        
        // Compute position and velocity of each particle j.
        // Euler's Method
        for(int j = 0; j < n; j++) {
            pstates[j].pos = pstates[j].vel * delta_t;
            pstates[j].vel = pstates[j].force * (delta_t / masses[j]);
        }
    }
    
    // outout
    for(int j = 0; j < n; j++) {
        cout << pstates[j].pos.x << " " << pstates[j].pos.y << " ";
        cout << pstates[j].vel.x << " " << pstates[j].vel.y << "\n";
    }
    return 0;
}
