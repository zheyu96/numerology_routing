#ifndef __PURIFICATION_H
#define __PURIFICATION_H

namespace Purification {

static inline double fidelity_to_werner(double fidelity) {
    return (4.0 * fidelity - 1.0) / 3.0;
}

static inline double werner_to_fidelity(double werner) {
    return (3.0 * werner + 1.0) / 4.0;
}

// Paper Eq. (6): P_p(w1, w2) = (1 + w1*w2) / 2.
static inline double success_prob_werner(double w1, double w2) {
    return (1.0 + w1 * w2) / 2.0;
}

// Paper Eq. (7): w_p(w1, w2) = (4*w1*w2 + w1 + w2) / (3*(1 + w1*w2)).
static inline double purified_werner(double w1, double w2) {
    return (4.0 * w1 * w2 + w1 + w2) / (3.0 * (1.0 + w1 * w2));
}

static inline double pumping_werner(double fresh_werner, int rounds) {
    double current_werner = fresh_werner;
    for(int r = 0; r < rounds; ++r) {
        current_werner = purified_werner(current_werner, fresh_werner);
    }
    return current_werner;
}

static inline double pumping_fidelity(double fresh_fidelity, int rounds) {
    double fresh_werner = fidelity_to_werner(fresh_fidelity);
    return werner_to_fidelity(pumping_werner(fresh_werner, rounds));
}

// Paper Eq. (8): Pr^(r+1) = Pr^(r) * Pr(u,v) * P_p(w^(r), w_e(u,v)).
static inline double pumping_success_prob(double entangle_prob, double fresh_werner, int rounds) {
    double current_werner = fresh_werner;
    double current_prob = entangle_prob;

    for(int r = 0; r < rounds; ++r) {
        current_prob *= entangle_prob;
        current_prob *= success_prob_werner(current_werner, fresh_werner);
        current_werner = purified_werner(current_werner, fresh_werner);
    }

    return current_prob;
}

} // namespace Purification

#endif
