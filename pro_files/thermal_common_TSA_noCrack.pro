Include "CvAg.pro";
Include "Ic.pro";
Include "thermal.data";

DefineConstant[
  N_ele = {1, Min 0, Max 20, Name "Parameters/# thin shell elements"}
];


//------------------------------------------------------------------------------
//----------------------------   Regions  --------------------------------------
//------------------------------------------------------------------------------
Group{
  dom_th_withoutDefect = Region [{}];
  dom_th = Region [{}];
  local_defect = Region[ {} ];

  For i In {1:4}
    tapes~{i} = Region[(2000 + i - 1)];
    bare_minus~{i} = Region[(20 + i - 1)];

    Shell += Region[{tapes~{i}}];
    dom_th_withoutDefect += Region[(20 + i - 1)];

    // add (potentially local defect)
    local_defect += Region[(40 + i - 1)];
    bare_minus~{i} += Region[(40 + i - 1)];
    dom_th += Region[{bare_minus~{i}}];
  EndFor

  bare_minus_5 = Region[bare_minus_1];

  bare_plus_1 = Region[bare_minus_3];
  bare_plus_2 = Region[bare_minus_4];
  bare_plus_3 = Region[bare_minus_1];
  bare_plus_4 = Region[bare_minus_2];
  bare_plus_5 = Region[bare_plus_1];

  tapes_0 = Region[tapes_4];
  tapes_5 = Region[tapes_1];

  /* For i In {1:4}
    dom_oneSide_minus~{i} = ElementsOf[ {tapes_helper~{i}, bare_helper~{i},
      bare_helper~{i+1}}, OnOneSideOf {tapes_helper~{i}}];

    dom_oneSide_plus~{i} = ElementsOf[ {tapes_helper~{i}, bare_helper_plus~{(i)},
      bare_helper_plus~{(i+1)}}, OnOneSideOf {tapes_helper~{i}}];
  EndFor */
}

  //------------------------------------------------------------------------------
  //---------------------------  Functions  --------------------------------------
  //------------------------------------------------------------------------------
  Function{

    w_Ins = w_Ins/(N_ele);

    // FIXME: assumed same layer thickness everywhere (okay since hom Neumann
    // for single layer ins anyways)
    kappa_stiffness[ ]
      = corr_fact * kKaptonStiffness[$1, $2]{order_1D_gauss, 2 * w_Ins};

    kappa_mass[ ] =corr_fact *
       kKaptonMass[$1, $2, $4, $5]{order_1D_gauss, 2 * w_Ins};

    heatCap_mass[ ] = corr_fact * CvKaptonMass[$1, $2, $4, $5]{order_1D_gauss, 2 * w_Ins};

    // kappa(T, B)
    kappa[dom_th] = (
                      f_hastelloy * kSteel_GetDP[$1] + f_silver * kAg_GetDP[$1] +
                      f_copper * kCu_NISTGetDP[$1,
                        rhoCu_NISTGetDP[$1, 0]{RRR_Cu, Tup_RRR_Cu},
                        rhoCu_NISTGetDP[$1, Norm[$2]]{RRR_Cu, Tup_RRR_Cu}]
                        {RRR_Cu} );

    // Cv(T)
    heatCap[dom_th] =  (
                          f_hastelloy * CvSteel_GETDP[$1] + f_silver * cvAg[$1] +
                          f_copper * CvCu_CUDIGetDP[$1]
                          );

    // rho
    rho_Copper[] = rhoCu_NISTGetDP[$1, Norm[$2]]{RRR_Cu, Tup_RRR_Cu};
    rho_Silver[] = rhoAG_CUDIGetDP[$1, Norm[$2]]{RRR_Cu, Tup_RRR_Cu};
    rho_Hastelloy = 1.25e-6;

    J[] = I/(cable_width * h_cable);
    Jcrit[] = Ic[$1]/(cable_width * h_cable);
	HTS_width = 1E-6;

    lambda[] = HomogenizedCC_CS_GetDP[Norm[$1], Norm[$2], $3, n[]] // J_cc, J_crit, rhoNC
      {Ecrit, absTol_CS, maxNumIter_CS, cable_width, HTS_width};

    rhoNC[] = 1./(f_hastelloy * 1./rho_Hastelloy +
        f_copper * 1./rho_Copper[$1,$2] +
        f_silver * 1./rho_Silver[$1, $2]);

    rhoHTS[] = Ecrit/$2 * (lambda[$1, $2, $3] * $1/$2)^(n-1);

    equivalentRho[] = $2 > 0? // critical current non-zero?
		(cable_width + HTS_width)/(cable_width/$3 + HTS_width/rhoHTS[$1, $2, $3]): // partial flow through SC
		(cable_width + HTS_width)/cable_width * $3; // no flow through SC
  }

//------------------------------------------------------------------------------
//----------------------   Jacobian and integration   --------------------------
//------------------------------------------------------------------------------

Jacobian{
  { Name jac_vol; // volume Jacobian
   Case{
  		{ Region All; Jacobian Vol; }
   }
  }

  { Name jac_sur ; // surface Jacobian
    Case {
        { Region All ; Jacobian Sur ;}
    }
  }
}

Integration{
	{ Name Int ; // Gauss integraion scheme
		Case {
			{ Type Gauss ;
				Case {
          { GeoElement Point;        NumberOfPoints  1; }
          { GeoElement Line;         NumberOfPoints  2; }
          { GeoElement Triangle;     NumberOfPoints  3; }
          { GeoElement Triangle2;     NumberOfPoints  3; }
          { GeoElement Triangle3;     NumberOfPoints  3; }

          { GeoElement Quadrangle;   NumberOfPoints  4; }
          { GeoElement Quadrangle2;   NumberOfPoints  4; }
          { GeoElement Tetrahedron;  NumberOfPoints  5; }
          { GeoElement Tetrahedron2;  NumberOfPoints  5; }
          { GeoElement Tetrahedron3;  NumberOfPoints  5; }

          { GeoElement Hexahedron;   NumberOfPoints  6; }
          { GeoElement Hexahedron2;   NumberOfPoints  6; }
          { GeoElement Hexahedron3;   NumberOfPoints  6; }

          { GeoElement Prism;        NumberOfPoints  6; }
          { GeoElement Pyramid;      NumberOfPoints  8; }
          { GeoElement Pyramid2;      NumberOfPoints  8; }

          { GeoElement Prism2;        NumberOfPoints  6; }

             }
			}
		}
	}
}
