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

#define M_SEED 123
#define M_ITER 10000



using namespace std;
using namespace std::chrono;

template<class T> void print_vector(vector<T>& vec)
{
    for(int i = 0; i < vec.size(); i++){
        cerr << vec[i] << " ";
    }
    cerr << endl;
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

public:
    Spaceship(int id = -1, int cs = -1, int np = -1, int nnp = -1, int fu = -1):
        GeneralSpaceVehicle(id, cs, np, nnp), m_followed_ufo(fu) {}

    int get_followed_ufo() { return m_followed_ufo; }

    void set_followed_ufo(int fu) { m_followed_ufo = fu; }

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
    int m_stars;
    int m_spaceships;
    int m_ufos;

    vector<Star> m_star_vector;
    vector<Spaceship> m_spaceship_vector;
    vector<Ufo> m_ufo_vector;

    default_random_engine m_engine;

public:

    StarTraveller();


    int init(vector<int> stars);
    vector<int> makeMoves(vector<int> ufos, vector<int> ships);

    void fill_spaceship_vector(vector<int> &ships);
    void update_spaceship_vector(vector<int> &ships);

    void fill_ufo_vector(vector<int> &ufos);
    void update_ufo_vector(vector<int> &ufos);


    void metropolis_ufo_spaceship_matching(vector<int> &destinations);
    void run_metropolis(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector);
    double get_distance_between_stars(Star &s1, Star &s2);
    double get_distance_between_spaceship_and_next_ufo_star(int spaceship_id, int ufo_id);
    double metropolis_energy(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector);
    void shuffle_parallel_vectors(uniform_int_distribution<int> &dist,vector<int> &ufos_parallel_vector);
    void shuffle_parallel_vectors(vector<int> &ufos_parallel_vector);

    void nearest_neighbour_ufo_spaceship_matching(vector<int> &destinations_vector);
    int find_next_nearest_ufo(int spaceship_id);

};


StarTraveller::StarTraveller(){
    m_stars = -1;
    m_spaceships = -1;
    m_ufos = -1;

    m_star_vector = vector<Star>();
    m_spaceship_vector = vector<Spaceship>();
    m_ufo_vector = vector<Ufo>();

    //srand(M_SEED);
    random_device rd;
    m_engine.seed(rd());
}


int StarTraveller::init(vector<int> stars) {

    m_stars = stars.size()/2;
    m_star_vector.resize(m_stars, Star());

    for(int i = 0; i < m_stars; i++){
        int x = stars[2*i];
        int y = stars[2*i + 1];
        m_star_vector[i].set_id(i);
        m_star_vector[i].set_x(x);
        m_star_vector[i].set_y(y);
    }

    return 0;
}


vector<int> StarTraveller::makeMoves(vector<int> ufos, vector<int> ships) {

    update_spaceship_vector(ships);
    update_ufo_vector(ufos);


    vector<int> destinations_vector(m_spaceships, 0);

    metropolis_ufo_spaceship_matching(destinations_vector);


    return destinations_vector;
}


void StarTraveller::fill_spaceship_vector(vector<int> &ships) {

    m_spaceships = ships.size();
    m_spaceship_vector.resize(m_spaceships, Spaceship());
    for(int i = 0; i < m_spaceships; i++)
        m_spaceship_vector[i].set_current_star( ships[i] );

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
    for(int i = 0; i < m_ufos; i++){
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
        for(int i = 0; i < m_ufos; i++){
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

    for(int i = 0; i < m_spaceships; i++)
        spaceships_parallel_vector[i] = i;

    for(int i = 0; i < m_ufos; i++)
        ufos_parallel_vector[i] = i;

    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    run_metropolis(spaceships_parallel_vector, ufos_parallel_vector);
    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //int elapsed_time = duration_cast<microseconds>( t2 - t1 ).count();
    //cerr << "Metropolis alg duration: " << elapsed_time << " ms." << endl;

    for(int i = 0; i < len; i++){
        if((spaceships_parallel_vector[i] != -1) && (ufos_parallel_vector[i] != -1)){

            destinations[ spaceships_parallel_vector[i] ] = m_ufo_vector[ ufos_parallel_vector[i] ].get_next_star();
            m_ufo_vector[ ufos_parallel_vector[i] ].set_following_spaceship(spaceships_parallel_vector[i]);
            m_spaceship_vector[ spaceships_parallel_vector[i] ].set_followed_ufo(ufos_parallel_vector[i]);

        } else if ((spaceships_parallel_vector[i] != -1) && (ufos_parallel_vector[i] == -1)) {
            destinations[ spaceships_parallel_vector[i] ] = m_spaceship_vector[ spaceships_parallel_vector[i] ].get_current_star();
        }

    }

}


void StarTraveller::run_metropolis(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector){

    uniform_int_distribution<int> dist(0, ufos_parallel_vector.size()-1);
    uniform_real_distribution<double> uni(0.0, 1.0);

    double d = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);
    vector<int> best_parallel_vector = ufos_parallel_vector;

    for(int i = 0; i < M_ITER; i++){
        vector<int> next_parallel_state = ufos_parallel_vector;
        //shuffle_parallel_vectors(next_parallel_state);
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
    cerr << "Optimized energy: " << d << endl;
    ufos_parallel_vector = best_parallel_vector;

}


double StarTraveller::get_distance_between_stars(Star &s1, Star &s2){

    double xd = s1.get_x() - s2.get_x();
    double yd = s1.get_y() - s2.get_y();
    return sqrt( xd*xd + yd*yd );

}


double StarTraveller::get_distance_between_spaceship_and_next_ufo_star(int spaceship_id, int ufo_id) {

    return get_distance_between_stars(m_star_vector[m_spaceship_vector[spaceship_id].get_current_star() ],
                                      m_star_vector[m_ufo_vector[ufo_id].get_next_star() ]);
}



double StarTraveller::metropolis_energy(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector){
    double energy = 0.0;
    for(int i = 0; i < spaceships_parallel_vector.size(); i++){
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


void StarTraveller::shuffle_parallel_vectors(vector<int> &ufos_parallel_vector) {

    int index_one = 0 + ( rand() % ( (ufos_parallel_vector.size()-1) - 0 + 1 ) );
    int index_two = 0 + ( rand() % ( (ufos_parallel_vector.size()-1) - 0 + 1 ) );

    int old_ufo_at_one = ufos_parallel_vector[index_one];

    ufos_parallel_vector[index_one] = ufos_parallel_vector[index_two];
    ufos_parallel_vector[index_two] = old_ufo_at_one;

}


void StarTraveller::nearest_neighbour_ufo_spaceship_matching(vector<int> &destinations_vector) {

    double travelled_distance = 0.0;
    for(int i = 0; i < m_spaceships; i++){

        int nearest_ufo = find_next_nearest_ufo(i);
        if(nearest_ufo != -1) {
            destinations_vector[i] = m_ufo_vector[nearest_ufo].get_next_star();
            m_ufo_vector[nearest_ufo].set_following_spaceship(i);
            m_spaceship_vector[i].set_followed_ufo(nearest_ufo);
            travelled_distance = travelled_distance +
                get_distance_between_stars(m_star_vector[m_spaceship_vector[i].get_current_star()], m_star_vector[m_ufo_vector[nearest_ufo].get_next_star()]);

            continue;
        } else {
            destinations_vector[i] = m_spaceship_vector[i].get_current_star();
            continue;
        }
    }
    cerr << "Travelled_distance : " << travelled_distance  << endl;

}


int StarTraveller::find_next_nearest_ufo(int spaceship_id){

    int ufo_id = -1;
    double min_distance = std::numeric_limits<double>::max();
    for(int i = 0; i < m_ufos; i++){
        if(m_ufo_vector[i].get_following_spaceship() == -1){
            double d = get_distance_between_stars( m_star_vector[ spaceship_id ], m_star_vector[ m_ufo_vector[i].get_next_star() ] );
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


int main()
{
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
}


/*
void StarTraveller::run_metropolis_swaps(vector<int> &spaceships_parallel_vector, vector<int> &ufos_parallel_vector){

    uniform_int_distribution<int> dist(0, ufos_parallel_vector.size()-1);
    uniform_real_distribution<double> uni(0.0, 1.0);

    double d = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);
    vector<int> best_parallel_vector = ufos_parallel_vector;

    for(int i = 0; i < M_ITER; i++){
        //vector<int> next_parallel_state = ufos_parallel_vector;

        double E1 = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);

        int index_one = dist(m_engine);
        int index_two = dist(m_engine);

        int old_ufo_at_one = ufos_parallel_vector[index_one];

        ufos_parallel_vector[index_one] = ufos_parallel_vector[index_two];
        ufos_parallel_vector[index_two] = old_ufo_at_one;

        double E2 = metropolis_energy(spaceships_parallel_vector, ufos_parallel_vector);
        double T = 1.0;

        double A = exp( (E1 - E2)/T );

        double p = uni(m_engine);

        if (p < A)
            //ufos_parallel_vector = next_parallel_state;
            if( E2 < d){
                d = E2;
                best_parallel_vector = ufos_parallel_vector;
            }
        else {
            int old_ufo_at_one = ufos_parallel_vector[index_one];

            ufos_parallel_vector[index_one] = ufos_parallel_vector[index_two];
            ufos_parallel_vector[index_two] = old_ufo_at_one;
        }

    }
    cerr << "Optimized energy: " << d << endl;
    ufos_parallel_vector = best_parallel_vector;

}
*/



