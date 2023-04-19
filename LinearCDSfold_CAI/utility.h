
#include "energy_parameter.h" // energy_parameter stuff

// pairs: 0:NP 1:CG 2:GC 3:GU 4:UG 5:AU 6:UA 7:NN
//  Vienna: 0:N 1:A 2:C 3:G 4:U 5:u 6:g

#include "intl11.h"
#include "intl21.h"
#include "intl22.h"

#define BASE(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: (x=='u'?5:(x=='g'?6:0)))))))
#define reBASE(x) ((x==1? 'A' : (x==2? 'C' : (x==3? 'G' : (x==4?'U': (x==5?'u':(x==6?'g':'*')))))))

#define MAXLOOP 30

inline int MIN2(int a, int b) {if (a <= b)return a;else return b;}
inline int MAX2(int a, int b) {if (a >= b)return a;else return b;}

#define NUC_TO_PAIR(x,y) (x==1? (y==4?5:0) : (x==2? (y==3?1:0) : (x==3 ? (y==2?2:(y==4?3:0)) : (x==4 ? (y==3?4:(y==1?6:0)) : 0))))
// #define NUC_TO_PAIR(x,y) (x==1? (y==4?5:(y=5?5:0)) : (x==2? (y==3?1:(y==6?1:0)) : (x==3? (y==2?2:(y==4?3:(y==5?3:0))) : (x==4? (y==3?4:(y==1?6:(y==5?4:0))) : (x==5? (y==3?4:(y==1?6:(y==5?4:0))) : (x==6 ? (y==2?2:(y==4?3:(y==5?3:0))) : 0))))))

inline int change_e_base(int base){
	if(base==BASE('u'))
    	return BASE('U');
    if(base==BASE('g'))
    	return BASE('G');
    return base;	
}
                       //AU Au                        //CG Cg                          //GC GU Gu                         //UG UA Ug                                  //uG uA ug                              //gC gU gu
inline int v_score_hairpin(int i, int j, int base_i, int base_ni, int base_pj, int base_j) {
    int size = j-i-1;

    base_i=change_e_base(base_i);
    base_ni=change_e_base(base_ni);
    base_pj=change_e_base(base_pj);
    base_j=change_e_base(base_j);

    int type = NUC_TO_PAIR(base_i, base_j);

    int energy;

    if(size <= 30)
        energy = hairpin37[size];
    else
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));

    if(size < 3) return energy; /* should only be the case when folding alignments */

    energy += mismatchH37[type][base_ni][base_pj];

    return energy;
}

// multi_loop
inline int E_MLstem(int type) {
    int energy = 0;

    
    if(type > 2) {
        energy += TerminalAU37;
    }

    energy += ML_intern37;

    return energy;
}

inline int v_score_multi(int nuci, int nucj) {

    nuci=change_e_base(nuci);
    nucj=change_e_base(nucj);

	int tt = NUC_TO_PAIR(nucj, nuci); // lhuang: closing pair in multi: reversed
    /* int si1 = NUM_TO_NUC(nuci1); */
    /* int sj1 = NUM_TO_NUC(nucj_1); */

    return E_MLstem(tt) + ML_closing37;
}

inline int v_score_M1(int nuci, int nuck) {

    nuci=change_e_base(nuci);
    nuck=change_e_base(nuck);

    int tt = NUC_TO_PAIR(nuci, nuck);
    /* int sp1 = NUM_TO_NUC(nuci_1); */
    /* int sq1 = NUM_TO_NUC(nuck1); */

    return E_MLstem(tt);

}

// exterior_loop
inline int v_score_external_paired(int nuci, int nucj) {

    nuci=change_e_base(nuci);
    nucj=change_e_base(nucj);

    int type = NUC_TO_PAIR(nuci, nucj);
    /* int si1 = NUM_TO_NUC(nuci_1); */
    /* int sj1 = NUM_TO_NUC(nucj1); */
    int energy = 0;

    

    if(type > 2)
        energy += TerminalAU37;
  return energy;
}

inline int v_score_external_unpaired(int i, int j) {
    return 0;
}

inline int v_score_single(int i, int j, int p, int q,
                        int nuci, int nuci1, int nucj_1, int nucj,
                        int nucp_1, int nucp, int nucq, int nucq1){

    nuci=change_e_base(nuci);
    nuci1=change_e_base(nuci1);
    nucj_1=change_e_base(nucj_1);
    nucj=change_e_base(nucj);

    nucp=change_e_base(nucp);
    nucp_1=change_e_base(nucp_1);
    nucq1=change_e_base(nucq1);
    nucq=change_e_base(nucq);

    int type = NUC_TO_PAIR(nuci, nucj);
    int type_2 = NUC_TO_PAIR(nucq, nucp);
    int n1 = p-i-1;
    int n2 = j-q-1;
    int nl, ns, u, energy;
    energy = 0;

    if (n1>n2) { nl=n1; ns=n2;}
    else {nl=n2; ns=n1;}

    if (nl == 0)
        return stack37[type][type_2];  /* stack */

    if (ns==0) {                      /* bulge */
        energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));
    if (nl==1) energy += stack37[type][type_2];
    else {
      if (type>2) energy += TerminalAU37;
      if (type_2>2) energy += TerminalAU37;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return int11_37[type][type_2][nuci1][nucj_1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = int21_37[type][type_2][nuci1][nucq1][nucj_1];
        else
          energy = int21_37[type_2][type][nucq1][nuci1][nucp_1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1];}
      else if (nl==3){              /* 2x3 loop */
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1];
    }
  }
  return energy;
}