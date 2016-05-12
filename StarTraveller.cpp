#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <memory>
#include <cmath>
#include <limits>

using namespace std;



template<class T> void print_vector(vector<T>& vec)
{
    for(int i = 0; i < vec.size(); i++){
        cerr << vec[i] << endl;
    }
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
    int m_if_ufo_present_id; // -1 if ufo is located at the star. Otherwise m_if_ufo_present_id gives the ufo id.

public:

    Star(int id = -1, int x = 0, int y = 0, bool status = false, int up = -1): m_id(id), m_x(x), m_y(y), m_status(status), m_if_ufo_present_id(up) {}
    ~Star() {}

    int get_id() const { return m_id; }
    int get_x() const { return m_x; }
    int get_y() const { return m_y; }
    bool get_status() const { return m_status; }
    int get_if_ufo_present_id() const { return m_if_ufo_present_id; }
    static int get_number_of_visited_stars() { return m_visited_stars; }

    void set_id(int id) { m_id = id; }
    void set_x(int x) { m_x = x; }
    void set_y(int y) { m_y = y; }
    void set_status(bool status);
    void set_as_visited();
    void set_if_ufo_present_id(int ufo_presence) { m_if_ufo_present_id = ufo_presence; }
};

int Star::m_visited_stars = 0;

void Star::set_status(bool status) {

    if(status == true){
        m_visited_stars = m_visited_stars + 1;
        m_status = status;
    } else {
        m_status = status;
    }

}


void Star::set_as_visited() {

    if(m_status == false){
        m_visited_stars = m_visited_stars + 1;
        m_status = true;
    }

}

ostream & operator<<(ostream & os, const Star &s){
    os << "Star id: " << s.get_id() << " x: " << s.get_x() << " y: " << s.get_y() << " status: " << s.get_status() << " Ufo presence: " << s.get_if_ufo_present_id();
    return os;
}



class Ufo {
private:
    int m_id;
    int m_current_star;
    int m_next_star;
    int m_next_to_next_star;
    bool m_is_occupied; // is the ufo followed by a spaceship

public:
    Ufo(int id = -1, int cs = -1, int ns = -1, int ntns = -1, bool is = false):
        m_id(id), m_current_star(cs), m_next_star(ns), m_next_to_next_star(ntns), m_is_occupied(is) {}
    ~Ufo() {}

    int get_id() const { return m_id; }
    int get_current_star() const { return m_current_star; }
    int get_next_star() const { return m_next_star; }
    int get_next_to_next_star() const { return m_next_to_next_star; }
    bool get_occupation() { return m_is_occupied; }

    void set_id(int id) { m_id = id; }
    void set_current_star(int cs) { m_current_star = cs; }
    void set_next_star(int ns) { m_next_star = ns; }
    void set_next_to_next_star(int ntns) { m_next_to_next_star = ntns; }
    void set_occupation(bool occup) { m_is_occupied = occup; }

};

ostream & operator<<(ostream & os, const Ufo &u){
    os << "Ufo id: " << u.get_id() << "\n Current star: " << u.get_current_star() << "\n Next star: " << u.get_next_star() << "\n Next to next star: " << u.get_next_to_next_star();
    return os;
}


class StarTraveller {

public:

    int n_stars;
    int n_ufos;
    int n_ships;

    vector<Star> v_stars;
    vector<Ufo> v_ufos;

    StarTraveller();


    int init(vector<int> stars);
    vector<int> makeMoves(vector<int> ufos, vector<int> ships);

    int find_nearest_star(int &csp);
    int find_nearest_ufo(int &csp);

    double get_distance_between_stars(Star &s1, Star &s2);
    void fill_the_ufo_vector(vector<int> &ufo_int, vector<Ufo> &ufo);
    void update_ufo_vector(vector<int> &ufo_int, vector<Ufo> &ufo);


    void set_the_ship_positions_to_ufo_positions(vector<Ufo> &v_ufos, vector<int> &ships);

};


StarTraveller::StarTraveller(){
    n_stars = -1;
    n_ufos = -1;
    n_ships = -1;
    v_stars = vector<Star>();
    v_ufos = vector<Ufo>();
}


int StarTraveller::init(vector<int> stars) {

    //cerr << "Number of stars: " << stars.size() << endl;

    n_stars = stars.size()/2;
    v_stars.resize(n_stars, Star());

    for(int i = 0; i < n_stars; i++){
        int x = stars[2*i];
        int y = stars[2*i + 1];
        v_stars[i].set_id(i);
        v_stars[i].set_x(x);
        v_stars[i].set_y(y);
    }


    return 0;
}

void StarTraveller::set_the_ship_positions_to_ufo_positions(vector<Ufo> &v_ufos, vector<int> &ships, vector<int> &v_destinations){

    if(n_ships > n_ufos){


        for(int ship_id = 0; ship_id < n_ufos; ship_id++){
            int ship_star = ships[ship_id]; // the star at which the ship is currently located

            // Ships start to follow the Ufos.
            int nearest_ufo = find_nearest_ufo(ship_star);
            if((nearest_ufo != -1) && (v_ufos[nearest_ufo].get_occupation() == false) && (v_ufos[nearest_ufo].get_current_star() != ship_star) ) {

                int destination = v_ufos[nearest_ufo].get_next_star();
                v_destinations[ship_id] = destination;
                v_stars[destination].set_as_visited();
                v_ufos[nearest_ufo].set_occupation(true);
                continue;
            }
        }

        for(int i = n_ufos; i < n_ships; i++)
            v_destinations[ship_id] = ships[ship_id];


    } else {




    }

}

vector<int> StarTraveller::makeMoves(vector<int> ufos, vector<int> ships) {

    n_ships = ships.size();
    vector<int> v_destinations(n_ships, 0);


    if(n_ufos == -1){
        n_ufos = ufos.size()/3;
        fill_the_ufo_vector(ufos, v_ufos);
    } else {
        update_ufo_vector(ufos, v_ufos);
    }



    for(int ship_id = 0; ship_id < n_ships; ship_id++){


        int ship_star = ships[ship_id]; // the star at which the ship is currently located

        cerr << "Number of stars: " << n_stars << endl;
        cerr << "Visited stars: " << Star::get_number_of_visited_stars() << endl;
        cerr << "Remaining stars: " << n_stars - Star::get_number_of_visited_stars() << endl;

        // Ships start to follow the Ufos.
        int nearest_ufo = find_nearest_ufo(ship_star);
        if((nearest_ufo != -1) && (v_ufos[nearest_ufo].get_occupation() == false) && (v_ufos[nearest_ufo].get_current_star() != ship_star) ) {

            int destination = v_ufos[nearest_ufo].get_next_star();
            v_destinations[ship_id] = destination;
            v_stars[destination].set_as_visited();
            v_ufos[nearest_ufo].set_occupation(true);
            continue;
        }


        /*

        int upi = v_stars[ship_star].get_if_ufo_present_id();
        if( upi != -1 ) {

            int next_ufo_star_destination = v_ufos[upi].get_next_star();
            //int next_to_next_ufo_star_destination = v_ufos[upi].get_next_to_next_star();

            //if( v_stars[next_ufo_star_destination].get_status() == false ){
            //if( ship_id != 0 ){
                v_destinations[ship_id] = next_ufo_star_destination;
                if( v_stars[next_ufo_star_destination].get_status() == false)
                    v_stars[next_ufo_star_destination].set_status(true);
                continue;
            //}
            //if( v_stars[next_to_next_ufo_star_destination].get_status() == false ){
            //    cerr << "Here" << endl;
            //    v_destinations[ship_id] = next_to_next_ufo_star_destination;
            //    v_stars[next_to_next_ufo_star_destination].set_status(true);
            //    continue;
            //}

        }



        int ns = find_nearest_star(ship_star);
        //cerr << "Closest star: " << ns << endl;

        if(ns != -1){
            v_destinations[ship_id] = ns;
            if(v_stars[ns].get_status() == false)
                v_stars[ns].set_status(true);
        } else {
            v_destinations[ship_id] = ship_star;
        }

        */

    }

    return v_destinations;
}



int StarTraveller::find_nearest_star(int &csp) {

    // csp - current ship position

    //Star current_star = v_stars[csp];
    int nsid = -1; // nsid - nearest stat id
    double dtns = std::numeric_limits<double>::max(); // dtns - distance to nearest star

    for(int i = 0; i < n_stars; i++){

        if(i == csp)
            continue;

        double d = get_distance_between_stars(v_stars[csp], v_stars[i]);
        if( (d < dtns) && (v_stars[i].get_status() == false) ){
            dtns = d;
            nsid = i;
        }
    }
    //cerr << "Closest distance: " << dtns << endl;
    return nsid;
}



int StarTraveller::find_nearest_ufo(int &csp) {

    // csp - current ship position
    int nus = -1; // nearest ufo star
    double dtnu = std::numeric_limits<double>::max(); // dtnu - distance to nearest ufo

    for(int i = 0; i < n_ufos; i++){
        double d = get_distance_between_stars(v_stars[csp], v_stars[v_ufos[i].get_current_star()]);
        if(d < dtnu){
            dtnu = d;
            nus = v_ufos[i].get_id();
        }
    }

    return nus;
}



double StarTraveller::get_distance_between_stars(Star &s1, Star &s2){

    double xd = s1.get_x() - s2.get_x();
    double yd = s1.get_y() - s2.get_y();
    return sqrt( xd*xd + yd*yd );

}


void StarTraveller::fill_the_ufo_vector(vector<int> &ufo_int_vector, vector<Ufo> &ufo_vector) {

    ufo_vector.resize(n_ufos, Ufo());

    for(int i = 0; i < n_ufos; i++){

        v_stars[ufo_int_vector[3*i]].set_if_ufo_present_id(i);

        ufo_vector[i].set_id(i);
        ufo_vector[i].set_current_star(ufo_int_vector[3*i]);
        ufo_vector[i].set_next_star(ufo_int_vector[3*i+1]);
        ufo_vector[i].set_next_to_next_star(ufo_int_vector[3*i+2]);
    }

}


void StarTraveller::update_ufo_vector(vector<int> &ufo_int_vector, vector<Ufo> &ufo_vector) {

    for(int i = 0; i < n_ufos; i++){
        int cs = ufo_vector[i].get_current_star();

        v_stars[cs].set_if_ufo_present_id(-1);
        v_stars[ufo_int_vector[3*i]].set_if_ufo_present_id(i);

        ufo_vector[i].set_current_star(ufo_int_vector[3*i]);
        ufo_vector[i].set_next_star(ufo_int_vector[3*i+1]);
        ufo_vector[i].set_next_to_next_star(ufo_int_vector[3*i+2]);
    }
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

