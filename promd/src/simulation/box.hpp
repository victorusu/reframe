#ifndef BOX_HPP
#define BOX_HPP

#include "../main.hpp"

struct Box {

    RVec len;

    // number invBoxX;
    // number invBoxY;
    // number invBoxZ;

    // number halfBoxX;
    // number halfBoxY;
    // number halfBoxZ;

    // void init()
    // {
    //     this->invBoxX = 1.0/this->x;
    //     this->invBoxY = 1.0/this->y;
    //     this->invBoxZ = 1.0/this->z;
    // }

    VHR_ALWAYS_INLINE
    number pbc(number r, const number boxaxis)
    {
        const number halfBox = 0.5 * boxaxis;

        while (r >  halfBox) r -= boxaxis;
        while (r < -halfBox) r += boxaxis;

        return r;
    }

    VHR_ALWAYS_INLINE
    RVec pbc(RVec rvec)
    {
        RVec r = rvec;
        const RVec halfBox = 0.5 * len;
        const RVec mhalfBox = - halfBox;

        // fprintf(stderr, " halfBox:  %f  %f  %f\n", halfBox[0], halfBox[1], halfBox[2]);
        // fprintf(stderr, "mhalfBox: %f %f %f\n", mhalfBox[0], mhalfBox[1], mhalfBox[2]);
        // fprintf(stderr, "    rvec: %f %f %f\n", rvec[0], rvec[1], rvec[2]);
        // fprintf(stderr, "       r: %f %f %f\n", r[0], r[1], r[2]);

        int i;

        for(i = 0; i < DIM; i++) {
            while (r[i] >  halfBox[i]) {
                // fprintf(stderr, "r[%d]: %f - -\n", i, r[i]);
                r[i] -= len[i];
            }
            while (r[i] <= mhalfBox[i]) {
                // fprintf(stderr, "r[%d]: %f - +\n", i, r[i]);
                r[i] += len[i];
            }
        }

        return r;
    }

    VHR_ALWAYS_INLINE
    number volume()
    {
    	return len[0] * len[1] * len[2];
    }

    VHR_ALWAYS_INLINE
    void scaleBy(const number mu)
    {
        len *= mu;
        // len *= mu;
        // len *= mu;
    }

    VHR_ALWAYS_INLINE
    number latticeShift(const number r, const number boxaxis)
    {
        const number invBox = 1.0 / boxaxis;
        const number halfBox = 0.5 * boxaxis;

        // Lattice shift
        const number shift = (r - halfBox) * invBox;

        return r - (nint(shift) * boxaxis);
    }


    VHR_ALWAYS_INLINE
    RVec latticeShift(const RVec rvec)
    {
        RVec halfBox = 0.5 * len;

        // Crap solution to the 1.0 / len
        RVec invBox;
        invBox = 1.0;
        invBox /= len;

        RVec shift = (rvec - halfBox) * invBox;

        return rvec - (rvecNINT(shift) * len);
    }

};

#endif //BOX_HPP