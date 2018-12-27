#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

typedef long long Z;

const double Pi = 4.0 * atan(1.0);
const double HalfPi = 2.0 * atan(1.0);
const double DegToRadCoef = Pi / 180.0;
const double Infty = numeric_limits<double>::infinity();

template <typename T>
int len(const T& cont) {
    return (int)cont.size();
}

struct Vec3 {
	double x;
	double y;
	double z;

	Vec3() {}
	Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

	static Vec3 fromLatLong(double la, double lo) {
		la *= DegToRadCoef;
		lo *= DegToRadCoef;
		double x = cos(lo) * cos(la);
		double y = sin(lo) * cos(la);
		double z = sin(la);
		return Vec3(x, y, z);
	}
};

double dot(Vec3 a, Vec3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3 cross(Vec3 a, Vec3 b) {
	return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
double norm2(Vec3 v) {
	return v.x * v.x + v.y * v.y + v.z * v.z;
}
double norm(Vec3 v) {
	return sqrt(norm2(v));
}
double sphereDist(Vec3 a, Vec3 b) {
	return 6378000.0 * (HalfPi - atan(dot(a, b) / norm(cross(a, b))));
}

struct Point {
    int child;
    Vec3 pos;
    Z wt;
    Z leftWt;
    Z rightWt;
};

const Z FinalWtLimit = 10000000;
const Point Korvatunturi = {-1, Vec3::fromLatLong(68.073611, 29.315278), 0, 0, 0};

struct Path {
    vector<Point> points;
    double cost;

    void finalize() {
        if(len(points) < 2 || points.front().child != -1 || points.back().child != -1) throw 0;
        Z wt = 0;
        for(int i = 0; i < len(points); ++i) {
            if(i != 0 && i != len(points) - 1 && points[i].child == -1) throw 0;
            wt += points[i].wt;
            points[i].leftWt = wt;
        }
        wt = 0;
        for(int i = len(points) - 1; i >= 0; --i) {
            wt += points[i].wt;
            points[i].rightWt = wt;
        }
        if(wt > FinalWtLimit) throw 0;
        cost = 0.0;
        for(int i = 1; i < len(points); ++i) {
            cost += sphereDist(points[i - 1].pos, points[i].pos);
        }
    }
};

vector<Path> initialPaths(Z wtLimit) {
    vector<Point> points;
	ifstream fp("nicelist.txt");
    if(!fp.good()) throw 0;
	while(true) {
		string item;
		fp >> item;
		if(!fp.good()) {
			break;
		}
		for(char& c : item) {
			if(c == ';') c = ' ';
		}
		item.push_back('\n');
		stringstream ss(item);
		ss.exceptions(ss.failbit | ss.badbit | ss.eofbit);
        Point point;
		double la, lo;
		ss >> point.child >> la >> lo >> point.wt;
        if(point.child == -1) throw 0;
		point.pos = Vec3::fromLatLong(la, lo);
        points.push_back(point);
    }

    Path path;
    path.points.push_back(Korvatunturi);

    while(!points.empty()) {
        double bestDist = Infty;
        int bestPi = -1;
        for(int pi = 0; pi < len(points); ++pi) {
            double dist = sphereDist(path.points.back().pos, points[pi].pos);
            if(dist < bestDist) {
                bestDist = dist;
                bestPi = pi;
            }
        }
        if(bestPi == -1) throw 0;
        swap(points[bestPi], points.back());
        path.points.push_back(points.back());
        points.pop_back();
    }

    path.points.push_back(Korvatunturi);

    while(true) {
        bool progress = false;
        for(int i = 1; i < len(path.points); ++i) {
            for(int j = i + 2; j < len(path.points); ++j) {
                double costDiff = 0.0;
                costDiff -= sphereDist(path.points[i - 1].pos, path.points[i].pos);
                costDiff -= sphereDist(path.points[j - 1].pos, path.points[j].pos);
                costDiff += sphereDist(path.points[i - 1].pos, path.points[j - 1].pos);
                costDiff += sphereDist(path.points[i].pos, path.points[j].pos);
                if(costDiff < 0.0) {
                    reverse(path.points.begin() + i, path.points.begin() + j);
                    progress = true;
                }
            }
        }
        if(!progress) break;
    }

    vector<Path> paths;
    int i = 1;
    while(i != len(path.points) - 1) {
        int j = i;
        Z wt = 0;
        while(j != len(path.points) - 1 && (j == i || wt + path.points[j].wt <= wtLimit)) {
            wt += path.points[j].wt;
            ++j;
        }
        Path sub;
        sub.points.push_back(Korvatunturi);
        sub.points.insert(sub.points.end(), path.points.begin() + i, path.points.begin() + j);
        sub.points.push_back(Korvatunturi);
        sub.finalize();
        paths.push_back(move(sub));
        i = j;
    }

    return paths;
}

void writeOutput(const vector<Path>& srcPaths) {
    vector<Path> paths(srcPaths.begin(), srcPaths.end());
    double cost = 0.0;
    for(Path& path : paths) {
        path.finalize();
        cost += path.cost;
    }
    Z score = (Z)round(cost);
    stringstream fnamess;
    fnamess << "output_" << score;
    string fname = fnamess.str();
    cout << "Writing " << fname << "\n";
    ofstream fp;
    fp.exceptions(fp.failbit | fp.badbit | fp.eofbit);
    fp.open(fname);
    for(const Path& path : paths) {
        bool first = true;
        for(const Point& point : path.points) {
            if(point.child != -1) {
                if(first) {
                    first = false;
                } else {
                    fp << ';';
                }
                fp << point.child;
            }
        }
        if(!first) {
            fp << '\n';
        }
    }
    fp.close();
}

mt19937 rng;

int main() {
	unsigned int seed = random_device{}();
	cout << "Using random seed " << seed << "\n";
	rng.seed(seed);

    double initTemp = exp(uniform_real_distribution<double>(log(5000.0), log(50000.0))(rng));
    cout << "Initial temperature " << initTemp << "\n";

    Z initWtLimit = 0;
    initWtLimit = max(initWtLimit, uniform_int_distribution<Z>(0, FinalWtLimit)(rng));
    initWtLimit = max(initWtLimit, uniform_int_distribution<Z>(0, FinalWtLimit)(rng));
    cout << "Initial weight limit " << initWtLimit << "\n";

    Z endWtLimit = 2 * FinalWtLimit;
    endWtLimit = min(endWtLimit, uniform_int_distribution<Z>(FinalWtLimit, 2 * FinalWtLimit)(rng));
    endWtLimit = min(endWtLimit, uniform_int_distribution<Z>(FinalWtLimit, 2 * FinalWtLimit)(rng));
    cout << "End weight limit " << endWtLimit << " (capped to " << FinalWtLimit << ")\n";

    long long iterations = (long long)3e11;
    cout << "Running " << (double)iterations << " iterations\n";

    cout << "Constructing initial paths\n";

    vector<Path> paths = initialPaths(initWtLimit);

    double cost = 0.0;
    for(const Path& path : paths) {
        cost += path.cost;
    }

    cout << "Initially " << len(paths) << " paths, cost = " << cost << "\n";

    for(long long iter = 0; iter < iterations; ++iter) {
        double prog = (double)iter / (double)iterations;
        auto temp = [&]() {
            double t = (1.0 - (double)iter / (double)iterations);
            return t * t * initTemp;
        };
        Z wtLimit = min(initWtLimit + iter * (endWtLimit - initWtLimit) / iterations, FinalWtLimit);

        if(!(iter & 0xFFFFFF)) {
            cout << "cost = " << cost << ", progress = " << prog << ", temperature = " << temp() << ", weight limit " << wtLimit << "\n";
        }
        if(!(iter & 0xFFF)) {
            int emptyCount = 0;
            int i = 0;
            while(i < len(paths)) {
                if(len(paths[i].points) == 2) {
                    if(emptyCount >= 200) {
                        swap(paths[i], paths.back());
                        paths.pop_back();
                        continue;
                    }
                    ++emptyCount;
                }
                ++i;
            }
            while(emptyCount < 200) {
                Path path;
                path.points.push_back(Korvatunturi);
                path.points.push_back(Korvatunturi);
                path.finalize();
                paths.push_back(move(path));
                ++emptyCount;
            }
        }

        int type = uniform_int_distribution<int>(0, 3)(rng);
        if(type == 0) {
            int ap = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            int bp = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            if(ap == bp) continue;
            int ai = uniform_int_distribution<int>(1, len(paths[ap].points) - 1)(rng);
            int bi = uniform_int_distribution<int>(1, len(paths[bp].points) - 1)(rng);

            if(paths[ap].points[ai - 1].leftWt + paths[bp].points[bi].rightWt > wtLimit) continue;
            if(paths[bp].points[bi - 1].leftWt + paths[ap].points[ai].rightWt > wtLimit) continue;

            double costDiff = 0.0;
            costDiff -= sphereDist(paths[ap].points[ai - 1].pos, paths[ap].points[ai].pos);
            costDiff -= sphereDist(paths[bp].points[bi - 1].pos, paths[bp].points[bi].pos);
            costDiff += sphereDist(paths[ap].points[ai - 1].pos, paths[bp].points[bi].pos);
            costDiff += sphereDist(paths[bp].points[bi - 1].pos, paths[ap].points[ai].pos);

            if(costDiff <= 0.0 || uniform_real_distribution(0.0, 1.0)(rng) < exp(-costDiff / temp())) {
                cost += costDiff;

                int ae = len(paths[ap].points);
                int be = len(paths[bp].points);
                paths[ap].points.insert(paths[ap].points.end(), paths[bp].points.begin() + bi, paths[bp].points.begin() + be);
                paths[bp].points.insert(paths[bp].points.end(), paths[ap].points.begin() + ai, paths[ap].points.begin() + ae);
                paths[ap].points.erase(paths[ap].points.begin() + ai, paths[ap].points.begin() + ae);
                paths[bp].points.erase(paths[bp].points.begin() + bi, paths[bp].points.begin() + be);

                paths[ap].finalize();
                paths[bp].finalize();
            }
        }
        if(type == 1) {
            int p = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            int i = uniform_int_distribution<int>(1, len(paths[p].points) - 1)(rng);
            int j = uniform_int_distribution<int>(1, len(paths[p].points) - 1)(rng);
            if(i > j) swap(i, j);
            if(j - i < 2) continue;

            double costDiff = 0.0;
            costDiff -= sphereDist(paths[p].points[i - 1].pos, paths[p].points[i].pos);
            costDiff -= sphereDist(paths[p].points[j - 1].pos, paths[p].points[j].pos);
            costDiff += sphereDist(paths[p].points[i - 1].pos, paths[p].points[j - 1].pos);
            costDiff += sphereDist(paths[p].points[j].pos, paths[p].points[i].pos);

            if(costDiff <= 0.0 || uniform_real_distribution(0.0, 1.0)(rng) < exp(-costDiff / temp())) {
                cost += costDiff;

                reverse(paths[p].points.begin() + i, paths[p].points.begin() + j);

                paths[p].finalize();
            }
        }
        if(type == 2) {
            int ap = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            int bp = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            if(ap == bp) continue;
            int as = uniform_int_distribution<int>(1, len(paths[ap].points) - 1)(rng);
            int ae = uniform_int_distribution<int>(1, len(paths[ap].points) - 1)(rng);
            int bs = uniform_int_distribution<int>(1, len(paths[bp].points) - 1)(rng);
            int be = uniform_int_distribution<int>(1, len(paths[bp].points) - 1)(rng);

            if(as > ae) swap(as, ae);
            if(bs > be) swap(bs, be);

            if(ae == as) continue;
            if(be == bs) continue;

            const Point& as1 = paths[ap].points[as - 1];
            const Point& as2 = paths[ap].points[as];
            const Point& ae1 = paths[ap].points[ae - 1];
            const Point& ae2 = paths[ap].points[ae];
            const Point& bs1 = paths[bp].points[bs - 1];
            const Point& bs2 = paths[bp].points[bs];
            const Point& be1 = paths[bp].points[be - 1];
            const Point& be2 = paths[bp].points[be];

            Z aWt1 = as1.leftWt;
            Z aWt2 = ae1.leftWt - aWt1;
            Z aWt3 = ae2.rightWt;
            Z bWt1 = bs1.leftWt;
            Z bWt2 = be1.leftWt - bWt1;
            Z bWt3 = be2.rightWt;

            if(aWt1 + bWt2 + aWt3 > wtLimit) continue;
            if(bWt1 + aWt2 + bWt3 > wtLimit) continue;

            double costDiff = 0.0;
            costDiff -= sphereDist(as1.pos, as2.pos);
            costDiff -= sphereDist(ae1.pos, ae2.pos);
            costDiff -= sphereDist(bs1.pos, bs2.pos);
            costDiff -= sphereDist(be1.pos, be2.pos);
            costDiff += sphereDist(as1.pos, bs2.pos);
            costDiff += sphereDist(ae1.pos, be2.pos);
            costDiff += sphereDist(bs1.pos, as2.pos);
            costDiff += sphereDist(be1.pos, ae2.pos);

            if(costDiff <= 0.0 || uniform_real_distribution(0.0, 1.0)(rng) < exp(-costDiff / temp())) {
                cost += costDiff;

                paths[ap].points.insert(paths[ap].points.begin() + ae, paths[bp].points.begin() + bs, paths[bp].points.begin() + be);
                paths[bp].points.insert(paths[bp].points.begin() + be, paths[ap].points.begin() + as, paths[ap].points.begin() + ae);
                paths[ap].points.erase(paths[ap].points.begin() + as, paths[ap].points.begin() + ae);
                paths[bp].points.erase(paths[bp].points.begin() + bs, paths[bp].points.begin() + be);

                paths[ap].finalize();
                paths[bp].finalize();
            }
        }
        if(type == 3) {
            int ap = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            int bp = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            int cp = uniform_int_distribution<int>(0, len(paths) - 1)(rng);
            if(ap == bp || ap == cp || bp == cp) continue;
            int ai = uniform_int_distribution<int>(1, len(paths[ap].points) - 1)(rng);
            int bi = uniform_int_distribution<int>(1, len(paths[bp].points) - 1)(rng);
            int ci = uniform_int_distribution<int>(1, len(paths[cp].points) - 1)(rng);

            if(paths[ap].points[ai - 1].leftWt + paths[bp].points[bi].rightWt > wtLimit) continue;
            if(paths[bp].points[bi - 1].leftWt + paths[cp].points[ci].rightWt > wtLimit) continue;
            if(paths[cp].points[ci - 1].leftWt + paths[ap].points[ai].rightWt > wtLimit) continue;

            double costDiff = 0.0;
            costDiff -= sphereDist(paths[ap].points[ai - 1].pos, paths[ap].points[ai].pos);
            costDiff -= sphereDist(paths[bp].points[bi - 1].pos, paths[bp].points[bi].pos);
            costDiff -= sphereDist(paths[cp].points[ci - 1].pos, paths[cp].points[ci].pos);
            costDiff += sphereDist(paths[ap].points[ai - 1].pos, paths[bp].points[bi].pos);
            costDiff += sphereDist(paths[bp].points[bi - 1].pos, paths[cp].points[ci].pos);
            costDiff += sphereDist(paths[cp].points[ci - 1].pos, paths[ap].points[ai].pos);

            if(costDiff <= 0.0 || uniform_real_distribution(0.0, 1.0)(rng) < exp(-costDiff / temp())) {
                cost += costDiff;

                int ae = len(paths[ap].points);
                int be = len(paths[bp].points);
                int ce = len(paths[cp].points);
                paths[ap].points.insert(paths[ap].points.end(), paths[bp].points.begin() + bi, paths[bp].points.begin() + be);
                paths[bp].points.insert(paths[bp].points.end(), paths[cp].points.begin() + ci, paths[cp].points.begin() + ce);
                paths[cp].points.insert(paths[cp].points.end(), paths[ap].points.begin() + ai, paths[ap].points.begin() + ae);
                paths[ap].points.erase(paths[ap].points.begin() + ai, paths[ap].points.begin() + ae);
                paths[bp].points.erase(paths[bp].points.begin() + bi, paths[bp].points.begin() + be);
                paths[cp].points.erase(paths[cp].points.begin() + ci, paths[cp].points.begin() + ce);

                paths[ap].finalize();
                paths[bp].finalize();
                paths[cp].finalize();
            }
        }
    }

    cout << "cost = " << cost << "\n";

    writeOutput(paths);

    return 0;
}
