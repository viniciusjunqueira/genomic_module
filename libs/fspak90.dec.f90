!Version of FSPAK90 for DEC; the DEC compiler crashes and pointer
!  variables are passed to old F77 programs and therefore this 
!  version avoids such passing by using allocatable variables; There
!  is some memory penalty but only during the initialization stage.


module sparseop
  contains

subroutine fspak90(operation,ija,rs4,sol4,det,msglev,maxmem,rank,rs8,sol8)
! interface to fspak; most parameters above optional and memory
! allocation dynamic
! VERSION with allocatable storage rather than pointer; works on DEC
  use sparsem
  implicit none
  
  interface
    subroutine fspak(mode,n,ia,ja,a,sol,flag,io1,io2,memory,&
                 mem_needed,work,i1,i2,i3,i4,i5,i6,rank)
    use kinds
    integer :: mode,n,ia(1),ja(1), flag,io1,io2,memory,&
                 mem_needed,i1,i2,i3,i4,i5,i6,rank 
    real (r8) :: a(1),sol(1),work(1)                       
    end subroutine
  end interface  
  
  character (len=*)  :: operation
  type (sparse_ija):: ija
  real(r8),optional :: rs8(:),sol8(:),det
  real,optional :: rs4(:),sol4(:)
  real (r8), allocatable:: sol(:)
 !optional parameters
  integer,optional::rank,msglev,maxmem
  integer :: msglev_o,maxmem_o
  integer,save :: rank_o,status=-1,last_n=0
        ! status = -1 (initial), 1 (ordering done), 2 (factorization done)
  integer,save :: memory=100        !initial memory, modified by later programs
  real (r8),save,allocatable :: work(:)
  real,external :: second

  ! for DEC copy the ija structure to allocatable
  integer::ijan,n1,n2
  integer,allocatable::ijaia(:),ijaja(:)
  real (r8),allocatable::ijaa(:)

  n1=size(ija%ia); n2=size(ija%ja)
  allocate(ijaia(n1),ijaja(n2),ijaa(n2))
  ijan=ija%n
  ijaia=ija%ia
  ijaja=ija%ja
  ijaa=ija%a



  if (last_n /= ija%n) then
        status=min(status,0)     !redo everything if dimension changes
        last_n=ija%n
  endif
!
  if (.not. present(msglev)) then
        msglev_o=0   ! default message level = minimal
      else
        msglev_o=msglev
  endif
  if (.not. present(maxmem)) then
        maxmem_o=huge(maxmem)         ! infinite memory limit if maxmem missing
     else
        maxmem_o=maxmem
  endif
  select case (operation(1:3))
       case ('fac')             !factorization
            status=min(status,1)        !nullify factorization if done
            call do_fact
       case('sol')            !solution
            allocate(sol(ija%n))
            if (present(sol8) .and. present(rs8)) then
                sol=rs8
                call do_fact
                call run_fspak('solve')
                sol8=sol
            elseif (present(sol4) .and. present(rs4)) then
                sol=rs4
                call do_fact
                call run_fspak('solve') 
                sol4=sol
              else
                print*,'FSPAK90: optional parameter SOL or RS missing'
                stop
            endif
            deallocate(sol)
       case('inv')           !sparse inversion
            call do_fact
            call run_fspak('inverse')
            status=1            !factorization lost
       case('res')            !reset storage
            if (status > -1) deallocate(work)
            status=-1
            memory=100
       case('det')              !determinant
            if (.not. present(det)) then
                print*,'FSPAK90: optional parameter DET missing'
              else
                call do_fact
                call run_fspak('det')
             endif
       case('lde')             !log determinant
            if (.not. present(det)) then
                print*,'FSPAK90: optional parameter DET misiing'
              else
                call do_fact
                call run_fspak('ldet')
            endif
       case('sta')             !print statistics
            call run_fspak('statistics')
       case default
          print*,'FSPAK90 operation ',operation,' unknown'
  end select

 ! copy back; only for DEC
  ija%ia=ijaia
  ija%ja=ijaja
  ija%a=ijaa
  deallocate(ijaia,ijaja,ijaa)


  contains
 
 subroutine do_fact
 ! executes all stages of fspak so that the factorization is available

 if (status == -1) then                
    allocate(work(memory/2+1))  ! real(r8) needed as integer
    status=0
 endif

 if (status == 0) then
    call run_fspak('ordering')
    call run_fspak('symbolic_fact')
    status=1
 endif
 
 if (status == 1) then
    call run_fspak('numeric_fact')
    status=2
 endif
 end subroutine do_fact


 subroutine run_fspak(op)
 ! runs operation op of fspak, adding more memory if necessary
 character (len=*) :: op
 real(r8),allocatable,save::work1(:)
 real(r8) :: sc(2)   !scratch
 integer :: mem_needed,i,flag

! print*,op
 do
   select case(op)
      case ('ordering')
            open(99,file='fspak90.ord')
            call fspak(10,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
            close (99)
            if (msglev_o > 2) then
               print*,'FSPAK90: size=',ija%n,'  #nonzeroes=',ija%nel
            endif
            if (msglev_o > 1) print*,'ordering time =',second()
      case('symbolic_fact')
            call fspak(20,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
            if (msglev_o > 1) print*,'symbolic factorization time=',second()
      case('numeric_fact')
            call fspak(40,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            if (msglev_o > 1) print*,'numerical factorization time=',second()
      case('solve')
            call fspak(50,ijan,ijaia,ijaja,ijaa,sol,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
      case('det')
            call fspak(54,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            det=sc(1)
      case('ldet')
            call fspak(55,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            det=sc(1)
      case('inverse')
            call fspak(61,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
             if (msglev_o >1) print*,'inversion time=',second()
      case('statistics')
            call fspak(80,ijan,ijaia,ijaja,ijaa,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
   end select
   if (flag == 0) then
          exit
      elseif (flag >= 4000000) then
        print*,'FSPAK: Zero or negative pivot for row ROW',&
            'during numerical factorization in row:',flag-4000000
        stop
     ! if memory insufficient, assign more
      elseif (flag > 100) then
        if (mem_needed > maxmem_o) then
            print*,'FSPAK90: memory needed:',mem_needed
            print*,'         maxmimum memory allowed:',maxmem
          else
            if (msglev_o>0) print*,'memory increased from:',memory,' to:',&
                                                          mem_needed
            allocate(work1(memory/2+1))     ! allocate memory for scratch
            work1=work                        !save old data
            deallocate(work)                  ! expand memory
            allocate(work(mem_needed/2+1))    !
            work(1:memory/2+1)=work1          !retrieve old data
            deallocate(work1)                 !release scratch
            memory=mem_needed
        endif
      else
        print*,'fspak: unknown error flag: ',flag
        stop
   endif
 enddo
 if (present(rank)) rank=rank_o
 
 end subroutine

end subroutine fspak90
 

end module
