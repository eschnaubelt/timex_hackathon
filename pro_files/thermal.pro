Include "thermal_common_TSA_noCrack.pro";

//------------------------------------------------------------------------------
//----------------------------   Constraints  ----------------------------------
//------------------------------------------------------------------------------

Constraint{
  { Name init_temp ;
    Case {
      { Region dom_th ; Type Init ; Value init_temp[] ; }
    }
  }

  { Name init_temp_shell ;
    Case {
      { Region Shell ; Type Init ; Value init_temp[] ; }
    }
  }
}

//------------------------------------------------------------------------------
//----------------------------   Function spaces  ------------------------------
//------------------------------------------------------------------------------

FunctionSpace{

  { Name Hgrad_T_tTinShell; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef Tn ; Function BF_Node ; Support dom_th ;
        Entity NodesOf[All, Not Shell]; }

      For i In {1:4}
        { Name sn_minus~{i} ; NameOfCoef Tn_minus~{i} ; Function BF_Node ;
          Support /*dom_oneSide_minus~{i}*/
          Region[{tapes~{i}, tapes~{i+1}, bare_minus~{i}, bare_minus~{i+1}}];
            Entity NodesOf[tapes~{i}, Not tapes~{i-1}] ; }

        { Name sn_plus~{i} ; NameOfCoef Tn_plus~{i} ; Function BF_Node ;
          Support /*dom_oneSide_plus~{i}*/
          Region[{tapes~{i}, tapes~{i+1}, bare_plus~{i}, bare_plus~{i+1}}] ;
          Entity NodesOf[tapes~{i}, Not tapes~{i-1}] ; }
      EndFor
    }

    SubSpace {
      { Name Shell_down; NameOfBasisFunction {sn_plus~{1}, sn_plus~{2},
        sn_plus~{3}, sn_plus~{4}}; }

      { Name Shell_up; NameOfBasisFunction {sn_minus~{1}, sn_minus~{2},
        sn_minus~{3}, sn_minus~{4}}; }

      /* { Name Tplus; NameOfBasisFunction {sn_plus~{1}, sn_plus~{2},
        sn_plus~{3}, sn_plus~{4}}; }

      { Name Tminus; NameOfBasisFunction {sn_minus~{1}, sn_minus~{2},
        sn_minus~{3}, sn_minus~{4}}; }

      { Name Tcont; NameOfBasisFunction {sn}; } */
    }

    Constraint {
        { NameOfCoef Tn ; EntityType NodesOf ; NameOfConstraint init_temp ; }

        For i In {1:4}
          { NameOfCoef Tn_minus~{i} ; EntityType NodesOf ;
              NameOfConstraint init_temp_shell ; }
          { NameOfCoef Tn_plus~{i} ; EntityType NodesOf ;
              NameOfConstraint init_temp_shell ; }
        EndFor


    }
  }

  For i In {1:N_ele-1}
    { Name Hgrad_T_tTinShell~{i} ; Type Form0 ;
      BasisFunction {
        { Name sn ; NameOfCoef Tn~{i} ; Function BF_Node ;
        Support Shell ; Entity NodesOf[ All ] ; }
      }

      Constraint {
        { NameOfCoef Tn~{i} ; EntityType NodesOf ;
          NameOfConstraint init_temp_shell ; }
      }
    }
  EndFor
}

//------------------------------------------------------------------------------
//----------------------------   Formulation(s)  -------------------------------
//------------------------------------------------------------------------------
Formulation {

  { Name thermal_form; Type FemEquation;
    Quantity {
      { Name T; Type Local; NameOfSpace Hgrad_T_tTinShell; }

      { Name Ti~{0}; Type Local; NameOfSpace Hgrad_T_tTinShell[Shell_down]; }
      For i In {1:N_ele-1}
        { Name Ti~{i} ; Type Local ; NameOfSpace Hgrad_T_tTinShell~{i}; }
      EndFor
      { Name Ti~{N_ele}; Type Local; NameOfSpace Hgrad_T_tTinShell[Shell_up]; }
    }

    Equation {
      Integral { [ kappa[{T}, magnField] * Dof{d T} , {d T} ] ;
        In dom_th; Integration Int ; Jacobian jac_vol ; }

      Integral { DtDof[heatCap[ {T} ] * Dof{T}, {T} ];
        In dom_th; Integration Int; Jacobian jac_vol;  }

      Integral { [- (cable_width + HTS_width)/cable_width * rhoNC[{T}, magnField] * SquNorm[J[]], {T} ];
        In local_defect; Integration Int; Jacobian jac_vol;  }

      Integral { [ - equivalentRho[J[], Jcrit[{T}], rhoNC[{T}, magnField]] *
        SquNorm[J[]], {T} ];
        In dom_th_withoutDefect; Integration Int; Jacobian jac_vol;  }

       For i In {0:N_ele-1}
         For k In {1:2} // row of the 1D FE matrix
           For l In {1:2} // column of the 1D FE matrix
            Integral {
                 [ kappa_mass[{Ti~{i}}, {Ti~{i+1}}, magnField, k, l] *
                 Dof{d Ti~{i + k - 1}} , {d Ti~{i + l - 1}}];
                 In Shell; Integration Int; Jacobian jac_sur;
             }

             Integral {
               [((k == l)? 1: -1) *
                 kappa_stiffness[{Ti~{i}}, {Ti~{i+1}}, magnField] *
                 Dof{Ti~{i + k - 1}} , {Ti~{i + l - 1}} ];
               In Shell; Integration Int; Jacobian jac_sur;
             }

              Integral {
               DtDof[  heatCap_mass[{Ti~{i}}, {Ti~{i+1}},
               magnField, k, l] * Dof{Ti~{i + k - 1}} , {Ti~{i + l - 1}} ];
               In Shell; Integration Int; Jacobian jac_sur;
             }

             EndFor // l: column of the 1D FE matrix
           EndFor // k: row of the 1D FE matrix
         EndFor
    }
  }
}
//------------------------------------------------------------------------------
//----------------------------   Resolution  -----------------------------------
//------------------------------------------------------------------------------

Resolution{
  { Name thermal_resolution;
    System {
      { Name thermal; NameOfFormulation thermal_form;
        NameOfMesh "../thermal_coil_TSA/thermal_coil_TSA_noCrack.msh"; }
    }

    Operation {

        SetExtrapolationOrder[0];

        InitSolution[thermal];

        PostOperation[resetMaxTemp];

        PostOperation[PrintMaxTemp];

        TimeLoopAdaptive[t_init, t_end,
                    t_step, t_step * 1E-6, 1, "BDF_2",
                    List[Breakpoints], System { { thermal, nl_relTol_time, nl_absTol_time, LinfNorm  } }
                  ]
                  {
                    IterativeLoopN[NMaxIt, relaxFactor,
                      System { { thermal, nl_relTol, nl_absTol, Solution LinfNorm } } ]
                        {
                          GenerateJac thermal ; SolveJac thermal ;
                        }
                    }
                    {
                        PostOperation[PrintMaxTemp];
                        Print[{#1}, Format "Max temperature is %g"];

                          Test[#1 > tempStop_2[]] {
                            Break[];
                          }

						Test[$Breakpoint >= 0] {
                              PostOperation[thermal];
                            //  Break[];
                          }
                    }

        // SaveSolution[thermal];
        // PostOperation[thermal];
    }
  }
}

//------------------------------------------------------------------------------
//--------------------------   Post-Processing  --------------------------------
//------------------------------------------------------------------------------
PostProcessing{

   { Name thermalPostPro; NameOfFormulation thermal_form;
     NameOfSystem thermal;
     Quantity {
       { Name T; Value {
           Term {[ {T} ]; In dom_th; Jacobian jac_vol; }
         }
       }

       { Name gradT; Value {
           Term {[ {d T} ]; In dom_th; Jacobian jac_vol; }
         }
       }

       { Name Tmax; Value{Term{ Type Global; [#1]; In dom_th;
         Jacobian jac_vol;} } }

       { Name Avg_temp; Value { Integral{ Type Global;
         [ {T}/(#2) ]; In dom_th;
         Jacobian jac_vol; Integration Int; } } }

         { Name Area; Value { Integral{ Type Global;
           [ 1 ]; In dom_th;
           Jacobian jac_vol; Integration Int; } } }
      }
    }
  }

//------------------------------------------------------------------------------
//--------------------------   Post-Operation  ---------------------------------
//------------------------------------------------------------------------------
PostOperation{
  { Name DoNoth; NameOfPostProcessing thermalPostPro; Operation {} }

  { Name thermal; NameOfPostProcessing thermalPostPro;
    Operation {
      /* Print[T, OnElementsOf dom_th, File "res_TSA/T_pre_allStep.pos"]; */

      Print[T, OnElementsOf dom_th, LastTimeStepOnly 1, Format Gmsh,
        OverrideTimeStepValue 0, File "res_TSA/TlastStep.pos",
        SendToServer "No"];

        /* Print[T, OnElementsOf dom_th, Format Gmsh,
          File "res_TSA/TlastStep.pos"]; */
    }
  }

  { Name PrintMaxTemp;  NameOfPostProcessing thermalPostPro;
      Operation {
        // Get maximum in bare region
        Print[ T, OnElementsOf dom_th, StoreMaxInRegister 1, Format Table ,
          LastTimeStepOnly 1, File "res_TSA/dummy.txt",
          SendToServer "No"] ;

        Print[Tmax, OnRegion dom_th, Format TimeTable, File StrCat["res_TSA/Tmax_combined",
          ".txt"], LastTimeStepOnly 1, AppendToExistingFile 1];

     /*   Print[T, OnElementsOf dom_th, File "res_TSA/T_combined.pos",
        AppendToExistingFile 1, LastTimeStepOnly 1]; */

       /* Print[gradT, OnElementsOf dom_th, File "res_TSA/gradT_combined.pos",
        AppendToExistingFile 1, LastTimeStepOnly 1]; */

        Print[Avg_temp, OnGlobal, File "res_TSA/Avg_temp_combined.txt",
          Format TimeTable,
          AppendToExistingFile 1, LastTimeStepOnly 1];

          Print[T, OnLine { { -center_radius - 2 * w_cable/4, 0, 0.001/2 } { -center_radius - 2 * w_cable/4 - no_turns * w_cable, 0, 0.001/2 } } { no_turns },
          File Sprintf["res_TSA/T_online_left%g.txt", w_Ins], Format TimeTable,
          AppendToExistingFile 1, LastTimeStepOnly 1
           ];

           Print[T, OnLine { { -1/Sqrt[2] * (center_radius + 3 * w_cable/8),
             1/Sqrt[2] * (center_radius + 3 * w_cable/8), 0.001/2 }
             { -1/Sqrt[2] * (center_radius + 3 * w_cable/8 + no_turns * w_cable),
              1/Sqrt[2] * (center_radius + 3 * w_cable/8 + no_turns * w_cable),
              0.001/2 } }
             { no_turns }, File Sprintf["res_TSA/T_online_upleft%g.txt", w_Ins], Format TimeTable,
             AppendToExistingFile 1, LastTimeStepOnly 1
            ];
      }
  }

  { Name resetMaxTemp;  NameOfPostProcessing thermalPostPro;
      Operation {

        // Print[ Q, OnElementsOf dom_th, File "res_TSA/Q.pos"] ;

        Echo["", Format Table, File StrCat["res_TSA/Tmax_combined", ".txt"] ];
        Echo["", Format Table, File StrCat["res_TSA/T_combined", ".pos"] ];
        Echo["", Format Table, File StrCat["res_TSA/Avg_temp_combined", ".txt"] ];
        Echo["", Format Table, File StrCat["res_TSA/T_online_left", Sprintf["%g", w_Ins], ".txt"] ];
        Echo["", Format Table, File StrCat["res_TSA/T_online_upleft", Sprintf["%g", w_Ins], ".txt"] ];

        Print[ Area, OnRegion dom_th, StoreInRegister 2, File "res_TSA/dummy.txt", Format Table ,
          LastTimeStepOnly 1] ;

      }
  }
}

/* DefineConstant[
  C_ = {"-sol -pos -v2 -ksp_rtol 1.e-12 -ksp_type gmres -ksp_gmres_restart 300 -pc_type ilu -pc_factor_levels 3", Name "GetDP/9ComputeCommand", Visible 0}
]; */

/* DefineConstant[
  C_ = {"-sol -pos -v2 -ksp_rtol 1.e-12  -ksp_type gmres  -ksp_gmres_restart 100  -pc_type jacobi  -pc_factor_levels 3", Name "GetDP/9ComputeCommand", Visible 0}
]; */
