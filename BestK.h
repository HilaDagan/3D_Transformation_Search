#ifndef BESTK_H
#define BESTK_H

#include <vector>
#include "Vector3.h"

/**
 * A transformation object, consisted of 6 degrees of freedom:
 * 3 of rotation and 3 of translation.
 */
struct transformation {
    double score = 0;

    // The translation in each of the coordinates:
    double x = -1;
    double y = -1;
    double z = -1;

    // The rotations:
    Rotation3 angles{};

    transformation(double score, double x, double y, double z, Rotation3 angles) :
            score(score), x(x), y(y), z(z), angles(angles) {}
};

/**
 * The vector is sorted in any moment according to the score, in a descending order.
 */
class BestK : public std::vector<transformation *> {
public:
    BestK(unsigned int k, bool toDel) : k_(k), toDel_(toDel) {};

    BestK() : k_(5), toDel_(true) {};

    void setK(unsigned int k) {
        k_ = k;
    }

    double score(transformation &t_model) const {
        return t_model.score;
    }

    bool push(transformation *in) {
        if (size() < k_) {
            push_back((transformation *) in);
            asort(in);
            return true;

        } else {
            transformation *atop = back();
            if (score(*atop) < score(*in)) {
                if (toDel_) {
                    delete atop;
                }
                back() = in;
                asort(in);
                return true;
            }
        }
        if (toDel_) {
            delete in;
        }
        return false;
    }

    void asort(transformation *in) {
        if (size() <= 1) return;

        reverse_iterator rfirst(rbegin());
        reverse_iterator rlast(rend());
        while (rfirst + 1 != rlast) {
            if (score(**(rfirst + 1)) >= score(*in)) {
                break;
            }
            *rfirst = *(rfirst + 1);
            *(rfirst + 1) = in;
            rfirst++;
        }
    }

    ~BestK() {
        if (toDel_) {
            for (unsigned int j = 0; j < size(); j++) {
                if (back() != nullptr) {
                    delete back();
                }
                pop_back();
            }
        }
    }

private:
    unsigned int k_; // the number of results we want
    bool toDel_;
};

#endif /* BESTK_H */
