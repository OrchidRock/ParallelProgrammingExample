#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <gmpxx.h>
#include <openmpi/mpi.h>

using namespace std;

#define ll long int
#define DIMENSION 2     // The dimension of Newton-Space.
#define MPF_STR_MAX_LEN 32
#define MAX_PARTICLE_NUMBER 100


#define GET_INDEX_PTR(p, i, j) (p + ((i)*DIMENSION + (j)) * MPF_STR_MAX_LEN)

const int BUFFER_MAX_LEN = MAX_PARTICLE_NUMBER * MPF_STR_MAX_LEN;
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
    
    friend ostream &operator<<(ostream &out,Vector &v) {
        out << v.x << " " << v.y;
        return out;
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


void get_input(int my_rank, char* masses_buf_ptr, char* pos_buf_ptr, 
               char* vel_buf_ptr, int* n_ptr, int *delta_ptr, int *T_ptr) {

    if(my_rank == 0) {
        int i, j, n, delta, T;
        cin >> n >> delta >> T; 
        *n_ptr = n;
        *delta_ptr = delta;
        *T_ptr = T;

        int local_n = min(n, MAX_PARTICLE_NUMBER);
        string tmp; 
        for(i = 0; i < local_n; i++) { 
            cin >> tmp;
            strncpy(masses_buf_ptr + i*MPF_STR_MAX_LEN, tmp.c_str(), 
                            min(int(tmp.size()), MPF_STR_MAX_LEN-1));
        }

        for(i = 0; i < local_n; i++) {
            for(j=0; j < DIMENSION; j++){
                cin >> tmp;
                strncpy(GET_INDEX_PTR(pos_buf_ptr, i, j), tmp.c_str(),
                            min(int(tmp.size()), MPF_STR_MAX_LEN-1));
            }
            for(j=0; j < DIMENSION; j++){
                cin >> tmp;
                strncpy(GET_INDEX_PTR(vel_buf_ptr, i, j), tmp.c_str(),
                            min(int(tmp.size()), MPF_STR_MAX_LEN-1));
            }
        }
    }
    MPI_Bcast(n_ptr, 1, MPI_INT, 0,  MPI_COMM_WORLD);
    MPI_Bcast(delta_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(T_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(masses_buf_ptr, BUFFER_MAX_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(pos_buf_ptr, BUFFER_MAX_LEN*DIMENSION, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(vel_buf_ptr, BUFFER_MAX_LEN*DIMENSION, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void pos_to_buffer(Vector &v, char* pos_buf_ptr, int index) {
    
    string x;
    string y;
    
    mp_exp_t exp;
    x = v.x.get_str(exp);
    x += "/";
    x += to_string(long(exp)); 
    
    y = v.y.get_str(exp);
    y += "/";
    y += to_string(long(exp));

    memset(GET_INDEX_PTR(pos_buf_ptr, index, 0), 0, MPF_STR_MAX_LEN);
    memset(GET_INDEX_PTR(pos_buf_ptr, index, 1), 0, MPF_STR_MAX_LEN);

    strncpy(GET_INDEX_PTR(pos_buf_ptr, index, 0), x.c_str(),
                    min(int(x.size()), MPF_STR_MAX_LEN-1));

    strncpy(GET_INDEX_PTR(pos_buf_ptr, index, 1), y.c_str(),
                    min(int(y.size()), MPF_STR_MAX_LEN-1));
}

/*
 * Reference: https://gmplib.org/list-archives/gmp-discuss/2002-December/000194.html.
 */
void buffer_to_pos(vector<Vector> &pos, char* pos_buf_ptr, int n) {
    int i;
    long long lExponent;
    long long pow_exp = 0;
    mpf_class tmp;
    for(i = 0; i < n; i++) {
        string x(GET_INDEX_PTR(pos_buf_ptr, i, 0));
        string y(GET_INDEX_PTR(pos_buf_ptr, i, 1));
        
        int elim_index = x.find('/');
        if(elim_index == string::npos) pos[i].x = x;
        else if (elim_index == 0) pos[i].x = 0;
        else{
            lExponent = stol(x.substr(elim_index+1));
            pos[i].x = x.substr(0, elim_index);
            pow_exp = lExponent - elim_index;
            pow_exp += (x[0] == '-') ? 1 : 0;
            tmp = pow(10.0, (double)pow_exp);
            pos[i].x *= tmp;
        }
        
        elim_index = y.find('/');
        if(elim_index == string::npos) pos[i].y = y;
        else if(elim_index == 0) pos[i].y = 0;
        else{
            lExponent = stol(y.substr(elim_index+1));
            pos[i].y = y.substr(0, elim_index);
            pow_exp = lExponent - elim_index;
            pow_exp += (y[0] == '-') ? 1 : 0;
            tmp = pow(10.0, (double)pow_exp);
            pos[i].y *= tmp;
        }
    }
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
    
    char vel_buffer[BUFFER_MAX_LEN*DIMENSION] = {'\0'};
    char pos_buffer[BUFFER_MAX_LEN*DIMENSION] = {'\0'};
    char masses_buffer[BUFFER_MAX_LEN] = {'\0'};
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

    get_input(my_rank, masses_buffer, pos_buffer, vel_buffer, &n, &delta, &T);
    
    vector<mpf_class> masses(n, 0); // the mass of particle.
    vector<Vector> pstates_pos(n);
    vector<Vector> pstates_vel(n);
    vector<Vector> pstates_force(n);
    
    for(int i=0; i<n; i++){
        masses[i] = masses_buffer+i*MPF_STR_MAX_LEN;
    }
    
    // The initial state
    buffer_to_pos(pstates_pos, pos_buffer, n);

    for(int i = 0; i < n; i++) {
        pstates_vel[i].x = GET_INDEX_PTR(vel_buffer, i, 0);
        pstates_vel[i].y = GET_INDEX_PTR(vel_buffer, i, 1);

        pstates_force[i]  = {0, 0};
    }
    
    mpf_class dist, dist_cubed;
    Vector diff(0,0);
    Vector force_jk(0,0);
    
    int local_n = n / comm_sz;
    int local_a = local_n*my_rank;
    int local_b = (my_rank != (comm_sz-1)) ? (local_a + local_n) : n;
    
     
    int local_buffer_len =  (local_n) * DIMENSION * MPF_STR_MAX_LEN;
    char *local_pos_buffer = new char [local_buffer_len]; 
    char *local_vel_buffer = new char [local_buffer_len]; 
    
    memset(local_pos_buffer, 0, local_buffer_len);
    memset(local_vel_buffer, 0, local_buffer_len);
     
    for(int i = 0; i < T; i++) {
        //cout << "i: " << i << endl;
        // Compute total force to each particle j.
        for(int j = local_a; j < local_b; j++) {
            pstates_force[j] = {0, 0};
        }
        
        for(int j = local_a; j < local_b; j++) {
            for(int k = 0; k < n; k++){
                
                if(k == j) continue;

                diff = pstates_pos[k] - pstates_pos[j]; 
                dist = diff.value();
                dist_cubed = dist*dist*dist;
                force_jk = diff * ((G * masses[j] * masses[k]) / dist_cubed);
                
                // Use Newton's  third law of motion.
                pstates_force[j] += force_jk;
                //pstates[k].force -= force_jk; 
            }
        }
        

        // Compute position and velocity of each particle j.
        // Euler's Method
        for(int j = local_a; j < local_b; j++) {
            pstates_pos[j] += (pstates_vel[j] * delta);
            pstates_vel[j] += pstates_force[j] * (delta / masses[j]); 

            int inx = j - local_a; 
            pos_to_buffer(pstates_pos[j], local_pos_buffer, inx);
        }

        MPI_Allgather(local_pos_buffer, local_buffer_len,
                   MPI_CHAR, pos_buffer, local_buffer_len, MPI_CHAR,MPI_COMM_WORLD);
        buffer_to_pos(pstates_pos, pos_buffer, n);
    }
    
    for(int j=local_a; j<local_b; j++) {
        pos_to_buffer(pstates_vel[j], local_vel_buffer, j-local_a);
    }

    MPI_Allgather(local_vel_buffer, local_buffer_len,
                    MPI_CHAR, vel_buffer, local_buffer_len, MPI_CHAR, MPI_COMM_WORLD);
    buffer_to_pos(pstates_vel, vel_buffer, n);

    // outout
    if(my_rank == 0) {
        cout << "Final State: " << endl; 
        for(int j = 0; j < n; j++) {
            cout << pstates_pos[j].x << " " << pstates_pos[j].y << " ";
            cout << pstates_vel[j].x << " " << pstates_vel[j].y << "\n";
        }
        cout << endl; 
    }
    delete [] local_pos_buffer;
    delete [] local_vel_buffer;
    MPI_Finalize();
    return 0;
}
