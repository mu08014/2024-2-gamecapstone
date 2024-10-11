using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Sorter
{
    private List<ulong> array_idx;
    private List<Tuple<int, int>> array_sup;

    private int ni;
    private int nj;
    private int nk;

    public Sorter(int ni_, int nj_, int nk_)
    {
        ni = ni_;
        nj = nj_;
        nk = nk_;
        resize(ni, nj, nk);
    }

    public void resize(int ni_, int nj_, int nk_)
    {
        array_sup = new List<Tuple<int, int>>(ni_ * nj_ * nk_);
        ni = ni_;
        nj = nj_;
        nk = nk_;
    }
}
