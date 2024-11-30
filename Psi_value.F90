Subroutine Psi_value(psi_flag,r_face,Courant,fi_small,psi_face)
    
    integer,intent(in) :: psi_flag
    Real,intent(in) :: r_face
    Real,intent(in) :: fi_small
    Real,intent(in) :: Courant
    Real :: fi_big
    Real,intent(out) :: psi_face

    If (psi_flag == 0) Then             !Upwind
        fi_big = 0.
        psi_face = fi_big
    ElseIf (psi_flag == 1) Then             !Lax-Wendroff
        fi_big = 1.d0
        psi_face = fi_big
    ElseIf (psi_flag == 2) Then             !Superbee
        fi_big = Max(fi_small,Max(Min(1.d0, 2.d0*r_face),Min(2.d0, r_face))) !Max(fi_small, Max(Min(1., 2.*r_face), Min(2., r_face)))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 3) Then         !Van Leer
        !fi_big = Max(fi_small, (r_face + abs(r_face)/(1 + abs(r_face)))) !HUGO
        fi_big = Max(fi_small, (r_face + abs(r_face))/(1 + r_face))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 4) Then         !MINMOD
        fi_big = Max(fi_small, Min(1., r_face))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 5) Then         !MUSCL
        fi_big = Max(fi_small, Min(2., Min(2.*r_face, (1. + r_face)/2.)))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 6) Then         !Super-C
        !
        If ((0. <= r_face) .and. (r_face <= 1.)) Then
            fi_big = Min(2*r_face/abs(Courant), 1.)
        ElseIf (r_face > 1.) Then
            fi_big = Min(r_face, 2./(1. - abs(Courant)))
        Else
            fi_big = fi_small
        EndIf
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 7) Then         !Ultimate-Quickest
        !Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
        fi_big = Max(fi_small, Max(0.5*(1 + r_face)+(1 - r_face)*(1 - 1*abs(Courant))/6., Max(2/(1 - abs(Courant)), 2*r_face/abs(Courant))))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 8) Then         !Hyper-C
        If (r_face > 0.) Then
            !Courant = 0.5*(HydroParam%Theta*HydroParam%u(iLayer,Face)*+(1.-HydroParam%Theta)*HydroParam%ut(iLayer,Face))*dt/MeshParam%CirDistance(Face)
            fi_big = Min(2*r_face/abs(Courant), 2/(1 - abs(Courant)))
        Else
            fi_big = fi_small
        EndIf
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 9) Then         !OSPRE
        fi_big = max(fi_small, 3*r_face*(r_face + 1)/(2*r_face**2 + 2*r_face + 2))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 10) Then        !SPL-1/3
        fi_big = max(fi_small, Min(2*r_face, Min(1./3. + 2*r_face/3., Min(2./3. + r_face/3., 2.))))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 11) Then        !H-CUI
        fi_big = max(fi_small, 3.*(r_face + abs(r_face))/2.*(r_face + 2))
        psi_face = fi_big - fi_small
    ElseIf (psi_flag == 12) Then        !SMART
        fi_big = max(fi_small, Min(2.*r_face, Min(3.*r_face/4. + 0.25, 4.)))
        psi_face = fi_big - fi_small
    EndIf    
    
    
    
    End Subroutine Psi_value
