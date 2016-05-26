#include <array>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <memory>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <deque>


#define M_ITER 30000
#define M_ITER_LONG 100000


#define M_ITER_LONG_1_100 32000000
#define M_ITER_LONG_100_200 18000000
#define M_ITER_LONG_200_300 13050000
#define M_ITER_LONG_300_400 9500000
#define M_ITER_LONG_400_500 7500000
#define M_ITER_LONG_500_600 6500000
#define M_ITER_LONG_600_700 5500000
#define M_ITER_LONG_700_800 4900000
#define M_ITER_LONG_800_900  3800000
#define M_ITER_LONG_900_1000 3100000
#define M_ITER_LONG_1000_INF 1500000



#define M_METROPOLIS_STAR_LIMIT 201

#define M_SPACESHIP_UFO_BALANCE 30
#define M_RADIOUS_MODIFIER 2.5

using namespace std;
using namespace std::chrono;

template<class T> void print_vector(vector<T>& vec)
{
    for(int i = 0; i < vec.size(); i++){
        cerr << vec[i] << " ";
    }
    cerr << endl;
}


template<class T> void print_vector(vector< vector<T> >& vec)
{
    for(int i = 0; i < vec.size(); i++){
        if(vec[i].size() == 0){
            cerr << "Empty";
        } else {
            for(int j = 0; j < vec[i].size(); j++){
                cerr << vec[i][j] << " ";
            }
        }
        cerr << endl;
    }
    cerr << endl;
}


template<class T> int number_of_elements_in_vector(vector< vector<T> >& vec)
{
    int N = 0;
    for(int i = 0; i < vec.size(); i++){
        N = N + vec[i].size();
    }
    return N;
}


template<class T> void getVector(vector<T>& v)
{
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}


class Star {
private:
    static int m_visited_stars;
    int m_id;
    int m_x;
    int m_y;
    bool m_status;
    int m_occupying_ufo; // -1 if no ufo is located at the star. Otherwise m_occupying_ufo gives the ufo id.

public:

    Star(int id = -1, int x = 0, int y = 0, bool status = false, int up = -1): m_id(id), m_x(x), m_y(y), m_status(status), m_occupying_ufo(up) {}
    ~Star() {}

    int get_id() const { return m_id; }
    int get_x() const { return m_x; }
    int get_y() const { return m_y; }
    bool get_status() const { return m_status; }
    int get_occupying_ufo() const { return m_occupying_ufo; }
    static int get_number_of_visited_stars() { return m_visited_stars; }

    void set_id(int id) { m_id = id; }
    void set_x(int x) { m_x = x; }
    void set_y(int y) { m_y = y; }
    void set_status(bool status) { m_status = status; }
    void set_as_visited();
    void set_occupying_ufo(int ufo_presence) { m_occupying_ufo = ufo_presence; }

};

int Star::m_visited_stars = 0;

void Star::set_as_visited() {

    if(m_status == false){
        m_visited_stars = m_visited_stars + 1;
        m_status = true;
    }

}


ostream & operator<<(ostream & os, const Star &s){
    os << "Star id: " << s.get_id() << " x: " << s.get_x() << " y: " << s.get_y() << " status: " << s.get_status() << " Ufo presence: " << s.get_occupying_ufo();
    return os;
}


class GeneralSpaceVehicle {
private:
    int m_id;
    int m_current_star;
    int m_next_star;
    int m_next_to_next_star;

public:
    GeneralSpaceVehicle(int id = -1, int cp = -1, int np = -1, int nnp = -1):
        m_id(id), m_current_star(cp), m_next_star(np), m_next_to_next_star(nnp) {}
    virtual ~GeneralSpaceVehicle() {}

    virtual int get_id() { return m_id; }
    virtual int get_current_star() { return m_current_star; }
    virtual int get_next_star() { return m_next_star; }
    virtual int get_next_to_next_star() { return m_next_to_next_star; }

    virtual void set_id(int id) { m_id = id; }
    virtual void set_current_star(int cp) { m_current_star = cp; }
    virtual void set_next_star(int np) { m_next_star = np; }
    virtual void set_next_to_next_star(int nnp) { m_next_to_next_star = nnp; }

};


class Spaceship : public GeneralSpaceVehicle {
private:
    int m_followed_ufo; // if -1 then the spaceship is not following any ufo
    int m_starting_star;

public:
    Spaceship(int id = -1, int cs = -1, int np = -1, int nnp = -1, int fu = -1, int ss = -1):
        GeneralSpaceVehicle(id, cs, np, nnp), m_followed_ufo(fu), m_starting_star(ss) {}

    int get_followed_ufo() { return m_followed_ufo; }
    int get_starting_star() { return m_starting_star; }

    void set_followed_ufo(int fu) { m_followed_ufo = fu; }
    void set_starting_star(int ss) { m_starting_star = ss; }

};


class Ufo : public GeneralSpaceVehicle {
private:
    int m_following_spaceship; // if -1 then the ufo is not followed

public:
    Ufo(int id = -1, int cs = -1, int np = -1, int nnp = -1, int fs = -1):
        GeneralSpaceVehicle(id, cs, np, nnp), m_following_spaceship(fs) {}

    int get_following_spaceship() { return m_following_spaceship; }

    void set_following_spaceship(int fs) { m_following_spaceship = fs; }
};


class StarTraveller {
private:
    int m_steps;
    int m_max_steps;
    int m_stars;
    int m_spaceships;
    int m_ufos;

    bool m_need_metropolis_ufo_matching;
    bool m_following_ufos;
    bool m_metropolis_with_reseting_execute_variable;
    bool m_need_combined_nearest_neigbour_and_metropolis_search;
    bool m_radious_wait;
    bool m_get_radii;

    int m_single_spaceship_counter;

    vector<int> m_spaceship_counter;

    vector<Star> m_star_vector;
    vector<Spaceship> m_spaceship_vector;
    vector<Ufo> m_ufo_vector;
    vector<double> m_radii_vector;

    default_random_engine m_engine;

    vector< vector<int> > m_many_spaceships_destination_vector;
    vector<int> m_single_destination_vector;

    vector<vector<double> > m_star_distances;

public:

    StarTraveller();


    int init(vector<int> stars);
    vector<int> makeMoves(vector<int> ufos, vector<int> ships);

    void fill_spaceship_vector(vector<int> &ships);
    void update_spaceship_vector(vector<int> &ships);

    void fill_ufo_vector(vector<int> &ufos);
    void update_ufo_vector(vector<int> &ufos);


    void metropolis_ufo_spaceship_matching(vector<int> &destinations);
    void run_metropolis_matching(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector);
    double get_distance_between_stars(int &s1, int &s2);
    double get_distance_between_stars(int &&s1, int &&s2);
    double get_distance_between_stars(int &s1, int &&s2);
    double get_distance_between_stars(int &&s1, int &s2);
    double get_distance_between_stars(Star &s1, Star &s2);
    double get_distance_between_spaceship_and_next_ufo_star(int spaceship_id, int ufo_id);
    double metropolis_energy(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector);
    void shuffle_parallel_vectors(uniform_int_distribution<int> &dist,vector<int> &ufos_parallel_vector);
    void shuffle_parallel_vectors(vector<int> &ufos_parallel_vector);

    void nearest_neighbour_ufo_spaceship_matching(vector<int> &destinations_vector);
    int find_next_nearest_ufo(int spaceship_id);

    void follow_ufo(vector<int> &destinations_vector);

    void move_to_nearest_star(vector<int> &destinations_vector);
    void make_nearest_neighbour_paths(bool add_spaceship);
    int get_closest_unvisited_star_to_spaceship(int &spaceship_id);
    double get_distance_to_closest_star_to_spaceship(int &spaceship_id);


    //void metropolis_get_destinations_for_every_spaceship_with_reseting(bool choice);
    void metropolis_get_destinations_for_every_spaceship_without_reseting(int number_of_iteration);
    //void metropolis_get_destinations_for_every_spaceship_single_ship();
    double get_full_spaceship_path_energy(int &spaceship_id, vector<int> &path);
    double get_full_spaceship_path_energy_omit_spaceship(int spaceship_id, vector<int> &path);
    double get_many_spaceships_destination_vector_energy();
    double get_many_spaceships_destination_vector_energy_omit_spaceship();
    double get_many_spaceships_destination_vector_energy(vector<int> &recalculate_vector, vector<int> &energy_vector);
    double get_many_spaceships_destination_vector_energy_omit_spaceship(vector<int> &recalculate_vector, vector<int> &energy_vector);


    void go_on_with_the_many_spaceships_moves(vector<int> &destinations_vector);
    void execute_metropolis_for_many_ships_with_reseting(int &&number_of_iteration, vector<int> &destinations_vector);
    void execute_combined_search(int &&number_of_iteration, vector<int> &destinations_vector);

    vector<int> two_opt_swap(vector<int> &v, int i, int j);
    void two_opt_optimize(int &spaceship_id, vector<int> &spaceship_path);
    void two_opt_optimize_on_all_paths();

    void radious_wait(vector<double> &radii, vector<int> &destinations_vector);
    bool check_if_all_ufos_are_followed();
    bool check_if_all_spaceships_are_following_ufos();
    vector<double> get_radii_for_each_spaceship();
    double distance_to_nearest_ufo(int &spaceship_id);

};


StarTraveller::StarTraveller(){
    m_steps = 0;
    m_max_steps = 0;
    m_stars = -1;
    m_spaceships = -1;
    m_ufos = -1;

    m_need_metropolis_ufo_matching = true;
    m_following_ufos = true;
    m_metropolis_with_reseting_execute_variable = true;
    m_need_combined_nearest_neigbour_and_metropolis_search = true;
    m_radious_wait = true;
    m_get_radii = true;

    m_single_spaceship_counter = 0;

    m_spaceship_counter = vector<int>();

    m_star_vector = vector<Star>();
    m_spaceship_vector = vector<Spaceship>();
    m_ufo_vector = vector<Ufo>();
    m_radii_vector = vector<double>();

    m_many_spaceships_destination_vector = vector< vector<int> >();
    m_single_destination_vector = vector<int>();

    m_star_distances = vector< vector<double> >();

    srand(time(0));

    random_device rd;
    m_engine.seed(rd());
}


int StarTraveller::init(vector<int> stars) {

    m_stars = stars.size()/2;
    m_star_vector.resize(m_stars, Star());

    m_max_steps = 4*m_stars;

    for(int i = 0; i < m_stars; i++){
        int x = stars[2*i];
        int y = stars[2*i + 1];
        m_star_vector[i].set_id(i);
        m_star_vector[i].set_x(x);
        m_star_vector[i].set_y(y);
    }

    m_star_distances.resize(m_stars);
    for(int i = 0; i < m_stars; ++i)
        m_star_distances[i].resize(m_stars);


    for(int i = 0; i < m_stars; ++i){
        for(int j = 0; j < i; ++j){
            double d = get_distance_between_stars( m_star_vector[i], m_star_vector[j]);
            m_star_distances[i][j] = d;
            m_star_distances[j][i] = d;
        }
    }

    return 0;
}


vector<int> StarTraveller::makeMoves(vector<int> ufos, vector<int> ships) {


    update_spaceship_vector(ships);
    update_ufo_vector(ufos);

    vector<int> destinations_vector(m_spaceships, 0);

    if(m_ufos == 0){
        execute_combined_search(0, destinations_vector);
        return destinations_vector;
    }

    int steps_left = m_max_steps - m_steps;
    int stars_left = m_stars - Star::get_number_of_visited_stars();

    if( (m_spaceships < M_SPACESHIP_UFO_BALANCE && m_ufos < M_SPACESHIP_UFO_BALANCE) ||
        (m_spaceships <= M_SPACESHIP_UFO_BALANCE && m_ufos > M_SPACESHIP_UFO_BALANCE) ){
        if(m_need_metropolis_ufo_matching == true) {

            cerr << "Ufo matching" << endl;
            metropolis_ufo_spaceship_matching(destinations_vector);
            m_need_metropolis_ufo_matching = false;

            m_steps++;
            return destinations_vector;
        } else {
            if( m_following_ufos == true) {

                follow_ufo(destinations_vector);
                m_steps++;
                if ( steps_left <= (stars_left+2) ){
                    cerr << endl;
                    cerr << "Max steps: " << m_max_steps << endl;
                    cerr << "Steps made: " << m_steps << endl;
                    cerr << "Steps left: " << steps_left << endl;
                    cerr << "Stars left: " << stars_left << endl;

                    m_following_ufos = false;
                }
                return destinations_vector;
            } else {
                execute_combined_search(M_ITER_LONG, destinations_vector);
                //cerr << "Number of steps in vector: " << number_of_elements_in_vector(m_many_spaceships_destination_vector) << endl;
            return destinations_vector;
            }
        }
    }


    if(m_radious_wait == true && steps_left > (stars_left+1)){

        if(m_get_radii == true){
            cerr << "Radious waiting" << endl;
            m_radii_vector = get_radii_for_each_spaceship();
            m_get_radii = false;
        }

        radious_wait(m_radii_vector, destinations_vector);
        m_steps++;
        return destinations_vector;
    } else {
        m_radious_wait = false;
        execute_combined_search(M_ITER_LONG, destinations_vector);
        m_steps++;
        return destinations_vector;
    }

}


void StarTraveller::fill_spaceship_vector(vector<int> &ships) {

    m_spaceships = ships.size();
    m_spaceship_vector.resize(m_spaceships, Spaceship());
    m_spaceship_counter = vector<int>(m_spaceships, 0);

    for(int i = 0; i < m_spaceships; i++){
        m_spaceship_vector[i].set_current_star( ships[i] );
        m_spaceship_vector[i].set_starting_star( ships[i] );
    }
}


void StarTraveller::update_spaceship_vector(vector<int> &ships){

    if(m_spaceships == -1){
        fill_spaceship_vector(ships);
    } else {
        for(int i = 0; i < m_spaceships; i++)
            m_spaceship_vector[i].set_current_star( ships[i] );
    }
}


void StarTraveller::fill_ufo_vector(vector<int> &ufos) {

    m_ufos = ufos.size()/3;
    m_ufo_vector.resize(m_ufos, Ufo());
    for(int i = 0; i < m_ufos; ++i){
        m_star_vector[ ufos[3*i] ].set_occupying_ufo(i);

        m_ufo_vector[i].set_id(i);
        m_ufo_vector[i].set_current_star( ufos[3*i] );
        m_ufo_vector[i].set_next_star( ufos[3*i+1] );
        m_ufo_vector[i].set_next_to_next_star( ufos[3*i+2] );
    }
}


void StarTraveller::update_ufo_vector(vector<int> &ufos) {

    if(m_ufos == -1) {
        fill_ufo_vector(ufos);
    } else {
        for(int i = 0; i < m_ufos; ++i){
            int cs = m_ufo_vector[i].get_current_star();

            m_star_vector[cs].set_occupying_ufo(-1);
            m_star_vector[ ufos[3*i] ].set_occupying_ufo(i);

            m_ufo_vector[i].set_id(i);
            m_ufo_vector[i].set_current_star( ufos[3*i] );
            m_ufo_vector[i].set_next_star( ufos[3*i+1] );
            m_ufo_vector[i].set_next_to_next_star( ufos[3*i+2] );
        }
    }
}


void StarTraveller::metropolis_ufo_spaceship_matching(vector<int> &destinations) {

    int len = (m_spaceships > m_ufos) ? m_spaceships : m_ufos;

    vector<int> spaceships_parallel_vector(len, -1);
    vector<int> ufos_parallel_vector(len, -1);

    for(int i = 0; i < m_spaceships; ++i)
        spaceships_parallel_vector[i] = i;

    for(int i = 0; i < m_ufos; ++i)
        ufos_parallel_vector[i] = i;

    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    run_metropolis_matching(spaceships_parallel_vector, ufos_parallel_vector);
    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //int elapsed_time = duration_cast<microseconds>( t2 - t1 ).count();
    //cerr << "Metropolis alg duration: " << elapsed_time << " ms." << endl;

    for(int i = 0; i < len; ++i){
        if((spaceships_parallel_vector[i] != -1) && (ufos_parallel_vector[i] != -1)){

            int next_ship_destination = m_ufo_vector[ ufos_parallel_vector[i] ].get_next_star();
            destinations[ spaceships_parallel_vector[i] ] = next_ship_destination;

            m_star_vector[next_ship_destination].set_as_visited();

            m_ufo_vector[ ufos_parallel_vector[i] ].set_following_spaceship(spaceships_parallel_vector[i]);
            m_spaceship_vector[ spaceships_parallel_vector[i] ].set_followed_ufo(ufos_parallel_vector[i]);

        } else if ((spaceships_parallel_vector[i] != -1) && (ufos_parallel_vector[i] == -1)) {
            destinations[ spaceships_parallel_vector[i] ] = m_spaceship_vector[ spaceships_parallel_vector[i] ].get_current_star();
        }

    }

}


void StarTraveller::run_metropolis_matching(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector){

    uniform_int_distribution<int> dist(0, ufos_parallel_vector.size()-1);
    uniform_real_distribution<double> uni(0.0, 1.0);

    double d = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);
    vector<int> best_parallel_vector = ufos_parallel_vector;

    for(int i = 0; i < M_ITER; ++i){
        vector<int> next_parallel_state = ufos_parallel_vector;
        shuffle_parallel_vectors(dist, next_parallel_state);

        double E1 = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);
        double E2 = metropolis_energy(spaceships_parallel_vector, next_parallel_state);
        double T = 1.0;

        double A = exp( (E1 - E2)/T );

        double p = uni(m_engine);

        if (p < A)
            ufos_parallel_vector = next_parallel_state;
            if( E2 < d){
                d = E2;
                best_parallel_vector = next_parallel_state;
            }
        else
            continue;
    }
    //cerr << "Optimized energy: " << d << endl;
    ufos_parallel_vector = best_parallel_vector;

}


inline double StarTraveller::get_distance_between_stars(int &s1, int &s2){

    return m_star_distances[s1][s2];
}


inline double StarTraveller::get_distance_between_stars(int &&s1, int &&s2){

    return m_star_distances[s1][s2];
}


inline double StarTraveller::get_distance_between_stars(int &s1, int &&s2){

    return m_star_distances[s1][s2];
}


inline double StarTraveller::get_distance_between_stars(int &&s1, int &s2){

    return m_star_distances[s1][s2];
}


double StarTraveller::get_distance_between_stars(Star &s1, Star &s2){

    double xd = s1.get_x() - s2.get_x();
    double yd = s1.get_y() - s2.get_y();
    return sqrt( xd*xd + yd*yd );

}


double StarTraveller::get_distance_between_spaceship_and_next_ufo_star(int spaceship_id, int ufo_id) {

    return m_star_distances[m_spaceship_vector[spaceship_id].get_current_star()][m_ufo_vector[ufo_id].get_next_star()];
}


double StarTraveller::metropolis_energy(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector){
    double energy = 0.0;
    for(int i = 0; i < spaceships_parallel_vector.size(); ++i){
        if((spaceships_parallel_vector[i] == -1) || (ufos_parallel_vector[i] == -1))
            continue;
        else
            energy = energy + get_distance_between_spaceship_and_next_ufo_star(spaceships_parallel_vector[i], ufos_parallel_vector[i]);
    }
    return energy;
}


void StarTraveller::shuffle_parallel_vectors(uniform_int_distribution<int> &dist, vector<int> &ufos_parallel_vector) {

    int index_one = dist(m_engine);
    int index_two = dist(m_engine);

    int old_ufo_at_one = ufos_parallel_vector[index_one];

    ufos_parallel_vector[index_one] = ufos_parallel_vector[index_two];
    ufos_parallel_vector[index_two] = old_ufo_at_one;

}


void StarTraveller::nearest_neighbour_ufo_spaceship_matching(vector<int> &destinations_vector) {

    double travelled_distance = 0.0;
    for(int i = 0; i < m_spaceships; ++i){

        int nearest_ufo = find_next_nearest_ufo(i);
        if(nearest_ufo != -1) {
            destinations_vector[i] = m_ufo_vector[nearest_ufo].get_next_star();
            m_ufo_vector[nearest_ufo].set_following_spaceship(i);
            m_spaceship_vector[i].set_followed_ufo(nearest_ufo);
            travelled_distance = travelled_distance +
                get_distance_between_stars(m_spaceship_vector[i].get_current_star(), m_ufo_vector[nearest_ufo].get_next_star());
            continue;
        } else {
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
            continue;
        }
    }
    //cerr << "Travelled_distance : " << travelled_distance  << endl;
}


int StarTraveller::find_next_nearest_ufo(int spaceship_id){

    int ufo_id = -1;
    double min_distance = std::numeric_limits<double>::max();
    for(int i = 0; i < m_ufos; ++i){
        if(m_ufo_vector[i].get_following_spaceship() == -1){
            double d = get_distance_between_stars( m_spaceship_vector[spaceship_id].get_current_star() ,  m_ufo_vector[i].get_next_star()  );
            if(d < min_distance){
                min_distance = d;
                ufo_id = i;
            }
        } else {
            continue;
        }
    }
    return ufo_id;
}


double StarTraveller::distance_to_nearest_ufo(int &spaceship_id){

    double min_distance = std::numeric_limits<double>::max();
    for(int i = 0; i < m_ufos; ++i){
            double d = get_distance_between_stars( m_spaceship_vector[spaceship_id].get_current_star() ,  m_ufo_vector[i].get_next_star()  );
            if(d < min_distance){
                min_distance = d;
            }
    }
    return min_distance;
}



void StarTraveller::follow_ufo(vector<int> &destinations_vector) {


    for(int i = 0; i < m_spaceships; ++i){
        int followed_ufo = m_spaceship_vector[i].get_followed_ufo();
        if(followed_ufo == -1){
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
        } else if (m_spaceship_vector[i].get_current_star() == m_ufo_vector[ followed_ufo ].get_next_to_next_star()) {
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
        } else {
            int next_ship_destination = m_ufo_vector[ followed_ufo ].get_next_star();
            destinations_vector[i] = next_ship_destination;
            m_star_vector[next_ship_destination].set_as_visited();
            m_spaceship_vector[i].set_starting_star(next_ship_destination);
        }
    }

}


void StarTraveller::move_to_nearest_star(vector<int> &destinations_vector) {

    for(int i = 0; i < m_spaceships; ++i){
        int closest_star = get_closest_unvisited_star_to_spaceship(i);
        if(closest_star != -1){
            destinations_vector[i] = closest_star;
            m_star_vector[closest_star].set_as_visited();
        } else {
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
        }
    }

}


void StarTraveller::make_nearest_neighbour_paths(bool add_spaceship = true) {

    m_many_spaceships_destination_vector.resize(m_spaceships, vector<int>());

    for(int i = 0; i < m_spaceships; ++i){
        int starting = m_spaceship_vector[i].get_starting_star();
        m_star_vector[starting].set_as_visited();
    }

    while(Star::get_number_of_visited_stars() < m_stars) {

        for(int i = 0; i < m_spaceships; ++i){
            int closest_star = get_closest_unvisited_star_to_spaceship(i);
            if(closest_star != -1){
                m_many_spaceships_destination_vector[i].push_back(closest_star);
                m_spaceship_vector[i].set_current_star(closest_star);
                m_star_vector[closest_star].set_as_visited();
            } else {
                if(add_spaceship)
                    m_many_spaceships_destination_vector[i].push_back(m_spaceship_vector[i].get_current_star());
            }
        }
    }

    for(int i = 0; i < m_spaceships; ++i) {
        int starting = m_spaceship_vector[i].get_starting_star();
        m_spaceship_vector[i].set_current_star(starting);
    }

}


int StarTraveller::get_closest_unvisited_star_to_spaceship(int &spaceship_id){

        int closes_star_id = -1;
        double closest_star_distance = std::numeric_limits<double>::max();
        int current_spaceship_star = m_spaceship_vector[ spaceship_id ].get_current_star();
        for(int j = 0; j < m_stars; ++j){
                if(m_star_vector[j].get_status() == false && j != m_spaceship_vector[spaceship_id].get_starting_star()) {
                    //double d = get_distance_between_stars( m_star_vector[current_spaceship_star] , m_star_vector[j] );
                    double d = get_distance_between_stars( current_spaceship_star , j );
                    if( d < closest_star_distance){
                        closest_star_distance  = d;
                        closes_star_id = j;
                    }
                }
        }
        return closes_star_id;

}


double StarTraveller::get_distance_to_closest_star_to_spaceship(int &spaceship_id){

        double closest_star_distance = std::numeric_limits<double>::max();
        int current_spaceship_star = m_spaceship_vector[ spaceship_id ].get_current_star();
        for(int j = 0; j < m_stars; ++j){
                if(m_star_vector[j].get_status() == false && j != m_spaceship_vector[spaceship_id].get_starting_star()) {
                    //double d = get_distance_between_stars( m_star_vector[current_spaceship_star] , m_star_vector[j] );
                    double d = get_distance_between_stars( current_spaceship_star , j );
                    if( d < closest_star_distance){
                        closest_star_distance  = d;
                    }
                }
        }
        return closest_star_distance;

}


void StarTraveller::metropolis_get_destinations_for_every_spaceship_without_reseting(int number_of_iteration) {

    // This function does not reset the m_many_spaceships_destination_vector
    // It uses its current state.
    // As a consequence it uses a somewhat different energy function since the spaceships
    // occupy stars that have not been visited.

    double energy = get_many_spaceships_destination_vector_energy_omit_spaceship();

    uniform_int_distribution<int> spaceships_dist(0, m_spaceships-1);
    uniform_real_distribution<double> uniform(0.0, 1.0);

    vector<vector<int> > best_paths = m_many_spaceships_destination_vector;
    double d = get_many_spaceships_destination_vector_energy_omit_spaceship();

    vector<int> recalculate_vector(m_spaceships, -1);
    vector<int> energy_vector(m_spaceships, 0.0);

    for(int i = 0; i < number_of_iteration; ++i){

        //if(i % 1000000 == 0)
        //   cerr << "We are at: " << i << endl;

        int spaceship_from = spaceships_dist(m_engine);
        int spaceship_to = spaceships_dist(m_engine);

        int n_stars_spaceship_from = m_many_spaceships_destination_vector[spaceship_from].size();
        int n_stars_spaceship_to = m_many_spaceships_destination_vector[spaceship_to].size();

        int from_offset = 0;
        int to_offset = 0;

        if(n_stars_spaceship_from == 0 && n_stars_spaceship_to == 0)
            continue;

        if(n_stars_spaceship_from == 0)
            continue;
        uniform_int_distribution<int> spaceship_from_star_dist(0, n_stars_spaceship_from-1);
        from_offset = spaceship_from_star_dist(m_engine);

        if(n_stars_spaceship_to == 0){
            to_offset = 0;
        } else {
            uniform_int_distribution<int> spaceship_to_star_dist(0, n_stars_spaceship_to-1);
            to_offset = spaceship_to_star_dist(m_engine);
        }

        double E1 = get_many_spaceships_destination_vector_energy_omit_spaceship(recalculate_vector, energy_vector);

        int spaceship_from_energy = energy_vector[spaceship_from];
        int spaceship_to_energy = energy_vector[spaceship_to];

        int transfered_star = 0;
        if(spaceship_from == spaceship_to){
            transfered_star = m_many_spaceships_destination_vector[spaceship_from][from_offset];
            m_many_spaceships_destination_vector[spaceship_from][from_offset] = m_many_spaceships_destination_vector[spaceship_to][to_offset];
            m_many_spaceships_destination_vector[spaceship_to][to_offset] = transfered_star;
        } else {
            auto pos_from = m_many_spaceships_destination_vector[spaceship_from].begin() + from_offset;
            auto pos_to = m_many_spaceships_destination_vector[spaceship_to].begin() + to_offset;

            transfered_star = (*pos_from);
            m_many_spaceships_destination_vector[spaceship_from].erase(pos_from);
            m_many_spaceships_destination_vector[spaceship_to].insert(pos_to, transfered_star);
        }

        recalculate_vector[spaceship_from] = -1;
        recalculate_vector[spaceship_to] = -1;

        double E2 = get_many_spaceships_destination_vector_energy_omit_spaceship(recalculate_vector, energy_vector);
        double T = 1.0;

        double A = exp( (E1 - E2)/T );

        double p = uniform(m_engine);

        if (p < A){
            if( E2 < d){
                d = E2;
                //cerr << "Recent minimum: " << d << endl;
                best_paths = m_many_spaceships_destination_vector;
            }
            continue;
        } else {

            if(spaceship_from == spaceship_to){
                transfered_star = m_many_spaceships_destination_vector[spaceship_from][from_offset];
                m_many_spaceships_destination_vector[spaceship_from][from_offset] = m_many_spaceships_destination_vector[spaceship_to][to_offset];
                m_many_spaceships_destination_vector[spaceship_to][to_offset] = transfered_star;

                energy_vector[spaceship_from] = spaceship_from_energy;
                energy_vector[spaceship_to] = spaceship_to_energy;

                recalculate_vector[spaceship_from] = 0;
                recalculate_vector[spaceship_to] = 0;

            } else {
                auto pos_from = m_many_spaceships_destination_vector[spaceship_from].begin() + from_offset;
                auto pos_to = m_many_spaceships_destination_vector[spaceship_to].begin() + to_offset;

                m_many_spaceships_destination_vector[spaceship_to].erase(pos_to);
                m_many_spaceships_destination_vector[spaceship_from].insert(pos_from, transfered_star);

                energy_vector[spaceship_from] = spaceship_from_energy;
                energy_vector[spaceship_to] = spaceship_to_energy;

                recalculate_vector[spaceship_from] = 0;
                recalculate_vector[spaceship_to] = 0;

            }
        }
    }

    m_many_spaceships_destination_vector = best_paths;

    //cerr << "After optimization" << endl;
    //cerr << "Dimension one: " << m_many_spaceships_destination_vector.size() << endl;
    energy = get_many_spaceships_destination_vector_energy_omit_spaceship();
    //cerr << "Full energy: " << energy << "\n\n";
}


double StarTraveller::get_full_spaceship_path_energy(int &spaceship_id, vector<int> &path){

    if(path.size() != 0){

        double path_energy = 0.0;
        path_energy = path_energy + get_distance_between_stars(m_spaceship_vector[spaceship_id].get_current_star(), path[0]);

        for(int i = 1; i <= path.size() - 1; ++i){
            path_energy = path_energy + get_distance_between_stars(path[i-1], path[i]);
        }
        return path_energy;
    } else {
        return 0.0;
    }

}


double StarTraveller::get_full_spaceship_path_energy_omit_spaceship(int spaceship_id, vector<int> &path){

    if(path.size() != 0){

        double path_energy = 0.0;
        path_energy = path_energy + get_distance_between_stars(m_spaceship_vector[spaceship_id].get_starting_star(), path[0]);

        for(int i = 1; i <= path.size() - 1; ++i){
            path_energy = path_energy + get_distance_between_stars(path[i-1], path[i]);
        }
        return path_energy;
    } else {
        return 0.0;
    }

}


double StarTraveller::get_many_spaceships_destination_vector_energy(){

    double energy = 0.0;
    for(int i = 0; i < m_many_spaceships_destination_vector.size(); ++i){
        cerr << "asd: " << get_full_spaceship_path_energy(i, m_many_spaceships_destination_vector[i]) << endl;
        energy = energy + get_full_spaceship_path_energy(i, m_many_spaceships_destination_vector[i]);
    }
    return energy;
}


double StarTraveller::get_many_spaceships_destination_vector_energy_omit_spaceship(){

    double energy = 0;
    for(int i = 0; i < m_many_spaceships_destination_vector.size(); ++i)
        energy = energy + get_full_spaceship_path_energy_omit_spaceship(i, m_many_spaceships_destination_vector[i]);
    return energy;
}


double StarTraveller::get_many_spaceships_destination_vector_energy(vector<int> &recalculate_vector, vector<int> &energy_vector) {

    double energy = 0.0;

    for(int i = 0; i < m_spaceships; i++){
        if(recalculate_vector[i] == -1){
            double en = get_full_spaceship_path_energy(i, m_many_spaceships_destination_vector[i]);
            energy_vector[i] = en;
            energy = energy + en;
            recalculate_vector[i] = 0;
        } else {
            energy = energy + energy_vector[i];
        }
    }

    return energy;
}


double StarTraveller::get_many_spaceships_destination_vector_energy_omit_spaceship(vector<int> &recalculate_vector, vector<int> &energy_vector) {

    double energy = 0;

    for(int i = 0; i < m_spaceships; i++){
        if(recalculate_vector[i] == -1){
            double en = get_full_spaceship_path_energy_omit_spaceship(i, m_many_spaceships_destination_vector[i]);
            energy_vector[i] = en;
            energy = energy + en;
        } else {
            energy = energy + energy_vector[i];
        }
    }

    return energy;
}


void StarTraveller::go_on_with_the_many_spaceships_moves(vector<int> &destinations_vector){

        for(int i = 0; i < m_spaceships; ++i){

            if(m_many_spaceships_destination_vector[i].size() != 0){
                if( m_spaceship_counter[i] < m_many_spaceships_destination_vector[i].size()) {
                    destinations_vector[i] = m_many_spaceships_destination_vector[i][ m_spaceship_counter[i] ];
                    m_spaceship_counter[i]++;
                } else {
                    destinations_vector[i] = m_spaceship_vector[i].get_current_star();
                }
            } else {
                destinations_vector[i] = m_spaceship_vector[i].get_current_star();
            }
        }
}


void StarTraveller::execute_metropolis_for_many_ships_with_reseting(int &&number_of_iteration, vector<int> &destinations_vector) {

    if(m_metropolis_with_reseting_execute_variable == true){
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        make_nearest_neighbour_paths(false);
        two_opt_optimize_on_all_paths();
        metropolis_get_destinations_for_every_spaceship_without_reseting(number_of_iteration);

        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        int elapsed_time = duration_cast<microseconds>( t2 - t1 ).count();
        cerr << "Metropolis alg duration: " << elapsed_time << " ms." << endl;

        m_metropolis_with_reseting_execute_variable = false;
    }
    go_on_with_the_many_spaceships_moves(destinations_vector);

}


void StarTraveller::execute_combined_search(int &&number_of_iteration, vector<int> &destinations_vector) {

    if(m_need_combined_nearest_neigbour_and_metropolis_search == true){

        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        make_nearest_neighbour_paths();
        double energy = get_many_spaceships_destination_vector_energy_omit_spaceship();
        cerr << "Get energy before two_opt: " << energy << endl;


        two_opt_optimize_on_all_paths();
        energy = get_many_spaceships_destination_vector_energy_omit_spaceship();
        cerr << "Get energy after two_opt: " << energy << endl;


        // We have different number of iteration for different number of stars so the
        // obtained answer is more accurate.
        if( m_stars >= 1 && m_stars < 100)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_1_100);
        if( m_stars >= 100 && m_stars < 200)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_100_200);
        if( m_stars >= 200 && m_stars < 300)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_200_300);
        if( m_stars >= 300 && m_stars < 400)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_300_400);
        if( m_stars >= 400 && m_stars < 500)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_400_500);
        if( m_stars >= 500 && m_stars < 600)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_500_600);
        if( m_stars >= 600 && m_stars < 700)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_600_700);
        if( m_stars >= 700 && m_stars < 800)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_700_800);
        if( m_stars >= 800 && m_stars < 900)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_800_900);
        if( m_stars >= 900 && m_stars < 1000)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_900_1000);
        if( m_stars >= 1000)
            metropolis_get_destinations_for_every_spaceship_without_reseting(M_ITER_LONG_1000_INF);


        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        int elapsed_time = duration_cast<microseconds>( t2 - t1 ).count();

        cerr << "Combined alg duration: " << elapsed_time << " ms." << endl;

        for(int i = 0; i < m_spaceships; ++i){
            destinations_vector[i] = m_spaceship_vector[i].get_starting_star();
        }

        m_need_combined_nearest_neigbour_and_metropolis_search = false;
    } else {
        go_on_with_the_many_spaceships_moves(destinations_vector);
    }
}


vector<int> StarTraveller::two_opt_swap(vector<int> &v, int i, int j) {

    vector<int> part_1(&v[0], &v[i]);
    vector<int> part_2(&v[i], &v[j+1]);
    reverse(part_2.begin(), part_2.end());

    vector<int> part_3(&v[j+1], &v[ v.size() ]);

    part_1.insert(part_1.end(), part_2.begin(), part_2.end());
    part_1.insert(part_1.end(), part_3.begin(), part_3.end());

    return part_1;
}


void StarTraveller::two_opt_optimize(int &spaceship_id, vector<int> &spaceship_path){

    int imp = 0;
    int path_size = spaceship_path.size();

    while(imp < 2){

        for(int i = 1; i < path_size-1; ++i){

            for(int j = i + 1; j < path_size; ++j){

                double d_before_prev = 0.0;
                double d_before_post = 0.0;

                double d_after_prev = 0.0;
                double d_after_post = 0.0;

                if(j != path_size - 1) {
                    d_before_prev = get_distance_between_stars( spaceship_path[i-1] , spaceship_path[i]);
                    d_before_post = get_distance_between_stars( spaceship_path[j] , spaceship_path[j+1]);

                    d_after_prev = get_distance_between_stars( spaceship_path[i-1] , spaceship_path[j] );
                    d_after_post = get_distance_between_stars( spaceship_path[i] , spaceship_path[j+1]);
                } else {
                    d_before_prev = get_distance_between_stars( spaceship_path[i-1] , spaceship_path[i]);
                    d_after_prev = get_distance_between_stars( spaceship_path[i-1] , spaceship_path[j] );
                }

                double d1 = d_before_prev + d_before_post;
                double d2 = d_after_prev + d_after_post;

                if( (d1-d2) > 0.0 ){
                    imp = 0;

                    vector<int> part_1(&spaceship_path[0], &spaceship_path[i]);
                    vector<int> part_2(&spaceship_path[i], &spaceship_path[j+1]);
                    reverse(part_2.begin(), part_2.end());
                    vector<int> part_3(&spaceship_path[j+1], &spaceship_path[ spaceship_path.size() ]);

                    part_1.insert(part_1.end(), part_2.begin(), part_2.end());
                    part_1.insert(part_1.end(), part_3.begin(), part_3.end());
                    spaceship_path = part_1;
                }
            }
        }
        imp++;

    }
}


void StarTraveller::two_opt_optimize_on_all_paths() {

    for(int i = 0; i < m_many_spaceships_destination_vector.size(); ++i)
        two_opt_optimize(i, m_many_spaceships_destination_vector[i]);
}


vector<double> StarTraveller::get_radii_for_each_spaceship(){

    vector<double> radii(m_spaceships, 0.0);
    vector<double> closest_ufos_vector(m_spaceships, 0.0);

    double avg = 0.0;
    for(int i = 0; i < m_spaceships; ++i){
        avg = avg + distance_to_nearest_ufo(i);
    }
    avg = avg/m_spaceships;

    for(int i = 0; i < m_spaceships; ++i){
        double closest_star_distance = get_distance_to_closest_star_to_spaceship(i);
        double closest_ufo_distance = distance_to_nearest_ufo(i);
        radii[i] = M_RADIOUS_MODIFIER*closest_star_distance;
    }
    return radii;
}


void StarTraveller::radious_wait(vector<double> &radii, vector<int> &destinations_vector){


    for(int i = 0; i < m_spaceships; ++i){
        double r = radii[i];
        if(m_spaceship_vector[i].get_followed_ufo() != -1){
            int ufo_id = m_spaceship_vector[i].get_followed_ufo();
            destinations_vector[i] = m_ufo_vector[ufo_id].get_next_star();
            m_star_vector[m_ufo_vector[ufo_id].get_next_star()].set_as_visited();
            continue;
        }

        vector<int> ufos_in_range;
        for(int j = 0; j < m_ufos; ++j){
            if(m_ufo_vector[j].get_following_spaceship() == -1){

                double d = get_distance_between_stars( m_spaceship_vector[i].get_current_star() ,  m_ufo_vector[j].get_next_star()  );
                if(d < r){
                    ufos_in_range.push_back(j);
                }
            } else {
                continue;
            }
        }

        double nearest_ufo_in_range = std::numeric_limits<double>::max();
        int nearest_ufo_id = -1;
        for(int k = 0; k < ufos_in_range.size(); ++k){
            int ufo_id = ufos_in_range[k];
            double d = get_distance_between_stars( m_spaceship_vector[i].get_current_star() ,  m_ufo_vector[ufo_id].get_next_star()  );

            if(d < nearest_ufo_in_range){
                nearest_ufo_in_range = d;
                nearest_ufo_id = ufo_id;
            }

        }

        if(nearest_ufo_id == -1){
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
            m_star_vector[m_spaceship_vector[i].get_current_star()].set_as_visited();
            continue;
        } else {
            destinations_vector[i] = m_ufo_vector[nearest_ufo_id].get_next_star();
            m_star_vector[m_ufo_vector[nearest_ufo_id].get_next_star()].set_as_visited();

            m_ufo_vector[nearest_ufo_id].set_following_spaceship(nearest_ufo_id);
            m_spaceship_vector[i].set_followed_ufo(nearest_ufo_id);
            continue;
        }
    }
}


bool StarTraveller::check_if_all_ufos_are_followed(){

    int followed_ufos = 0;
    for(int i = 0; i < m_ufos; i++){
        if(m_ufo_vector[i].get_following_spaceship() != -1)
            followed_ufos++;
    }

    if(followed_ufos == m_ufos)
        return true;
    else
        return false;

}


bool StarTraveller::check_if_all_spaceships_are_following_ufos(){

    int followed_spaceships = 0;
    for(int i = 0; i < m_spaceships; i++){
        if(m_spaceship_vector[i].get_followed_ufo() != -1)
            followed_spaceships++;
    }

    if(followed_spaceships == m_spaceships)
        return true;
    else
        return false;
}



int main()
{

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    int NStars;
    cin >> NStars;
    vector<int> stars(NStars);
    getVector(stars);

    StarTraveller algo;
    int ignore = algo.init(stars);
    cout << ignore << endl;
    cout.flush();

    while (true)
    {
        int NUfo;
        cin >> NUfo;
        if (NUfo<0) break;
        vector<int> ufos(NUfo);
        getVector(ufos);
        int NShips;
        cin >> NShips;
        vector<int> ships(NShips);
        getVector(ships);
        vector<int> ret = algo.makeMoves(ufos, ships);
        cout << ret.size() << endl;
        for (int i = 0; i < ret.size(); ++i) {
            cout << ret[i] << endl;
        }
        cout.flush();
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    int elapsed_time = duration_cast<microseconds>( t2 - t1 ).count();
    //cerr << "Execution time " << elapsed_time << " ms. which is " << 100.0*elapsed_time/20000000.0 << " %." << endl;

}

