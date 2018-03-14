//
//  qsort.h
//  stap++
//
//  Created by six on 2017/12/21.
//

#ifndef qsort_h
#define qsort_h

struct sortnode{
    int key;
    double value;
};


inline unsigned int Partition(unsigned int *a,unsigned int lo,unsigned int hi){
    unsigned int c=lo+rand()%(hi-lo+1);
    unsigned int pivot=a[c];
    a[c]=a[lo];
    a[lo]=pivot;
    while(lo<hi){
        while((lo<hi)&&(pivot<=a[hi]))hi--;
        a[lo]=a[hi];
        while((lo<hi)&&(a[lo]<=pivot))lo++;
        a[hi]=a[lo];
    }
    a[lo]=pivot;
    return lo;
}

void QuickSort(unsigned int *a,unsigned int lo, unsigned int hi){
    if(hi-lo<2)return;
    unsigned int mi=Partition(a,lo,hi-1);
    QuickSort(a,lo,mi);
    QuickSort(a,mi+1,hi);
}



inline int Partition(sortnode *a,int lo,int hi){
    int c=lo+rand()%(hi-lo+1);
    sortnode pivot=a[c];
    a[c]=a[lo];
    a[lo]=pivot;
    while(lo<hi){
        while((lo<hi)&&(pivot.key<=a[hi].key))hi--;
        a[lo]=a[hi];
        while((lo<hi)&&(a[lo].key<=pivot.key))lo++;
        a[hi]=a[lo];
    }
    a[lo]=pivot;
    return lo;
}

inline void QuickSort(sortnode* a,int lo, int hi){
    if(hi-lo<2)return;
    int mi=Partition(a,lo,hi-1);
    QuickSort(a,lo,mi);
    QuickSort(a,mi+1,hi);
}



inline long long Partition(long long *a,double* b,long long lo,long long hi){
    //a is key, while b is value
    long long c=lo+rand()%(hi-lo+1);
    long long pivot=a[c];
    double pivot_b=b[c];
    a[c]=a[lo];  b[c]=b[lo];
    a[lo]=pivot;  b[lo]=pivot_b;
    while(lo<hi){
        while((lo<hi)&&(pivot<=a[hi]))hi--;
        a[lo]=a[hi]; b[lo]=b[hi];
        while((lo<hi)&&(a[lo]<=pivot))lo++;
        a[hi]=a[lo]; b[hi]=b[lo];
    }
    a[lo]=pivot;
    b[lo]=pivot_b;
    return lo;
}

void QuickSort(long long* a,double* b,long long lo, long long hi){
    if(hi-lo<2)return;
    long long mi=Partition(a,b,lo,hi-1);
    QuickSort(a,b,lo,mi);
    QuickSort(a,b,mi+1,hi);
}


#endif /* qsort_h */
