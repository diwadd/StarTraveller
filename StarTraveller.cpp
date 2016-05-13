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



class Star {
private:
    int m_id;
    int m_x;
    int m_y;
    bool m_status;

public:

    Star(int id = -1, int x = 0, int y = 0, bool status = false): m_id(id), m_x(x), m_y(y), m_status(status) {}
    ~Star() {}

    int get_id() const { return m_id; }
    int get_x() const { return m_x; }
    int get_y() const { return m_y; }
    bool get_status() const { return m_status; }

    void set_id(int id) { m_id = id; }
    void set_x(int x) { m_x = x; }
    void set_y(int y) { m_y = y; }
    void set_status(bool status) { m_status = status; }
};


ostream & operator<<(ostream & os, const Star &s){
    os << "Star id: " << s.get_id() << " x: " << s.get_x() << " y: " << s.get_y() << " status: " << s.get_status();
    return os;
}



class UFO {
private:
    int m_current_star;
    int m_next_star;
    int m_next_to_next_star;

public:
    UFO(int cs = -1, int ns = -1, int ntns = -1): m_current_star(cs), m_next_star(ns), m_next_to_next_star(ntns) {}
    ~UFO() {}

    int get_current_star() { return m_current_star; }
    int get_next_star() { return m_next_star; }
    int get_next_to_next_star() { return m_next_to_next_star; }

    void set_current_star(int cs) { m_current_star = cs; }
    void set_next_star(int ns) { m_next_star = ns; }
    void set_next_to_next_star(int ntns) { m_next_to_next_star = ntns; }

};



class StarTraveller {

public:

    int NStars;
    vector<int> used;

    int n_stars;
    vector<Star> v_stars;
    //const static double g_max_double;


    int init(vector<int> stars)
    {
        NStars = stars.size()/2;
        used.resize(NStars, 0);


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

    vector<int> makeMoves(vector<int> ufos, vector<int> ships)
    {

        int n_ships = ships.size();
        vector<int> v_destinations(n_ships, 0);

        //cerr << "Number of ships " << n_ships << endl;
        //cerr << "destinations vector size " << v_destinations.size() << endl;

        for(int ship_id = 0; ship_id < n_ships; ship_id++){

            int ns = find_nearest_star(ships[ship_id]);
            //cerr << "Closest star: " << ns << endl;

            if(ns != -1){
                v_destinations[ship_id] = ns;
                v_stars[ns].set_status(true);
            } else {
                v_destinations[ship_id] = ships[ship_id];
            }

        }

        return v_destinations;


    }


    int find_nearest_star(int &csp);
    double get_distance_between_stars(Star &s1, Star &s2);

};



int StarTraveller::find_nearest_star(int &csp) {

    // csp - current ship position

    Star current_star = v_stars[csp];
    int nsid = -1; // nsid - nearest stat id
    double dtns = std::numeric_limits<double>::max(); // dtns - distance to nearest star

    for(int i = 0; i < n_stars; i++){

        if(i == csp)
            continue;

        double d = get_distance_between_stars(current_star, v_stars[i]);
        if( (d < dtns) && (v_stars[i].get_status() == false) ){
            dtns = d;
            nsid = i;
        }
    }
    //cerr << "Closest distance: " << dtns << endl;
    return nsid;
}


double StarTraveller::get_distance_between_stars(Star &s1, Star &s2){

    double xd = s1.get_x() - s2.get_x();
    double yd = s1.get_y() - s2.get_y();
    return sqrt( xd*xd + yd*yd );

}


// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v)
{
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}


//const double StarTraveller::g_max_double = std::numeric_limits<double>::max();

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

