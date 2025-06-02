
#include <R.h>

static int is_unrelated(int *ipa, int *ima, int i, int j)
{
    if (i == j)
        return 0;
    if (ipa[i] < 0 && ipa[j] < 0)
        return 1;
    if (ipa[i] < 0 || j < i) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    return is_unrelated(ipa, ima, ipa[i], j)
        & is_unrelated(ipa, ima, ima[i], j);
}

// pedigree has n individuals, numbered 0 to n - 1
// ipa gives the father of each individual, either one of the individuals
//     in the pedigree (hence integer between 0 and n - 1) or a negative
//     number to indicate missing.
// ima gives the mother of each individual coded the same way as fathers
// offspring must come before parents, that is i < ipa[i] and i < ima[i]
// result is either TRUE or FALSE

void aster_validtrips(int *ipa, int *ima, int *nindin, int *result)
{
    int nind = nindin[0];

    for (int i = 0; i < nind; ++i) {
        if (ipa[i] >= nind || ima[i] >= nind)
            error("ipa or ima out of range");
        if ((ipa[i] < 0) != (ima[i] < 0))
            error("every individual must have two parents or none");
        if ((ipa[i] >= i) != (ima[i] >= i))
            error("offspring must come before parents");
    }

    result[0] = 1;
    for (int i = 0; i < nind; ++i)
        if  (ipa[i] >= 0)
            result[0] &= is_unrelated(ipa, ima, ipa[i], ima[i]);
}

