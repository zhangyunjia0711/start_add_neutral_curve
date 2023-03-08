!------------------------------------------------------------------------------
! TJU/Department of Mechanics, Fluid Mechanics, Code START
!------------------------------------------------------------------------------
!> @Author: Jianxin Liu shookware@tju.edu.cn
!> @Date: 2022-09-05 09:32:33
!> @LastEditTime: 2023-03-08 11:48:29
!> @LastEditors: Jianxin Liu shookware@tju.edu.cn
!> @Description: 
!> @FilePath: /START_all/src/public/module_cmd_progress.f90
!> @Copyright (c) 2023 by Jianxin Liu email: shookware@tju.edu.cn, All Rights Reserved.
!------------------------------------------------------------------------------
Module CMD_Progress

    use stdlib_string_type
    Implicit None
    private
    Logical, parameter, public :: CMD_PROGRESS_ABSOLUTE = .true.
    Type, public :: CLS_CMD_Progress
        Integer, private :: N, lens, i
        Character :: M = "#", O = "."
        type(string_type) :: Prefix
    Contains
        Procedure :: Set
        Procedure :: Put
    End Type CLS_CMD_Progress

contains

    Subroutine Set(this, N, L)
        Class(CLS_CMD_Progress) :: this
        Integer, Intent(IN) :: N, L
        this%N = N
        this%lens = L
        this%i = 0
        this%Prefix = " Progress: "
    End Subroutine Set

    Subroutine Put(this, K, bAbsol)
        Class(CLS_CMD_Progress) :: this
        Integer, Intent(IN) :: K
        Logical, optional :: bAbsol
        Character(len=1) :: br
        integer :: jm
        this%i = this%i + K
        if (present(bAbsol)) then
            if (bAbsol) this%i = K
        end if
        if (this%i > this%n) this%i = this%n
        jm = Nint(real(this%i*this%lens)/real(this%N))
        if (this%i < this%n) then
            br = char(13)
        else
            br = char(10)
        end if
        write (*, '(5a,f6.2,2a)', advance="no") this%Prefix//'', '[', &
            ! write( * , '(5a,f6.2,2a\)') trim(this%Prefix) , '[' , & !// 如您的编译器不支持，请用上方语句代替
            repeat(this%M, jm), repeat(this%O, this%lens - jm), '] ', this%i*100.0/this%N, "%", br
    End Subroutine Put

End Module CMD_Progress
