#include <iostream>

#define N 100


struct tour_t {
    int cities[N];
    int n;
    int cost;    
};

struct my_stack_t {
    int list[(N*N)/2];
    int list_sz;
};

int Tour_city(tour_t tour, int i) {
    return tour.cities[i];
}

//#define Tour_city(tour, i) (tour.cities[i])

void Push(my_stack_t stack, int city) {
    int loc = stack.list_sz;
    stack.list[loc] = city;
    stack.list_sz++;
}

void Push_copy(my_stack_t stack, tour_t tour) {
    
}

int main(int argc, char *argv[])
{
        
    return 0;
}
