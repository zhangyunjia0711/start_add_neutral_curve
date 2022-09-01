module sort_functions

   use iso_fortran_env

   implicit none

   private

   integer, parameter :: rs_kind = real32, rd_kind = real64

  !! This module provides sorting routines for integer, real and
  !! character arrays. A generic interface is used, so the same
  !! commands will work with all three types of input.
  !!
  !! There is a functional interface (.sort.) which returns an
  !! an array which, if it is used to index the original array,
  !! will produce an array sorted in ascending order. There is
  !! also a subroutine interface, which will sort the array that
  !! is passed to it and return it. The subroutines just call
  !! the functional interface and use the result to index the
  !! array.
  !!
  !! There is also a convenience function to return an index of
  !! the unique elements in an array. This is slightly harder
  !! to use as straight subscripting will not work -- the index
  !! array will be zero padded with the number of zeroes depending
  !! on how many duplicate values there are.

  !! $Log: sort_functions.f90,v $
  !! Revision 1.3  2005/10/14 04:24:39  aidan
  !! Changed the basic sort 'engine' to a quicksort after unacceptably long
  !! sorting times for 2000x2000 synchrotron data.
  !!
  !! Revision 1.2  2003/09/12 03:52:23  aidan
  !! Added support for sorting double precision lists. Added a reverse function.
  !!
  !! Revision 1.1  2003/08/19 07:03:55  aidan
  !! Initial revision
  !!

   ! Revision control software updates this character parameter.
   ! The 'ident' command can extract this version string from an
   ! object file or executable, which means one can identify which
   ! version of the module was used to compile it.
   character(len=*), parameter :: version = "$Id: sort_functions.f90,v 1.3 2005/10/14 04:24:39 aidan Exp aidan $"

   ! A generic interface to the operator which returns sorted indices
   interface operator(.unique.)
      module procedure unique_integer_mask, unique_real_mask, unique_double_mask, unique_string_mask
   end interface

   ! A generic interface to the operator which returns sorted indices
   interface num_unique
      module procedure num_unique_integers, num_unique_reals, num_unique_doubles, num_unique_strings
   end interface

   ! A generic interface to the operator which returns sorted indices
   interface unique
      module procedure unique_integers, unique_reals, unique_doubles, unique_strings
   end interface

   ! A generic interface to the operator which returns sorted indices
   interface operator(.sort.)
      ! module procedure sort_real_return_key, sort_double_return_key, sort_integer_return_key, sort_strings_return_keys
      module procedure qsortd_integer, qsortd_real, qsortd_double, sort_strings_return_keys
   end interface

   ! A subroutine interface -- calls operator interface and then reorders
   ! input array
   interface sort
      module procedure sort_integer, sort_real, sort_double, sort_string
   end interface

   ! Simple function to reverse the order of a 1D array
   interface reverse
      module procedure reverse_integer, reverse_real, reverse_double, reverse_string
   end interface

   !
   ! Public stuff
   !

   ! Public routines
   public :: sort, unique, num_unique, reverse

   ! Public operators
   public :: operator(.sort.), operator(.unique.)

contains

   subroutine sort_integer(list)

      ! Performs a heapsort into ascending order of the input list -- calls
      ! another routine to provide a sorted index key, and then just sets it to this

      ! Input and output variable
      integer, dimension(:), intent(inout) :: list

      list = list(.sort.list)

      return
   end subroutine sort_integer

   subroutine sort_real(list)

      ! Performs a heapsort into ascending order of the input list -- calls
      ! another routine to provide a sorted index key, and then just sets it to this

      ! Input and output variable
      real, dimension(:), intent(inout) :: list

      list = list(.sort.list)

      return
   end subroutine sort_real

   subroutine sort_double(list)

      ! Performs a heapsort into ascending order of the input list -- calls
      ! another routine to provide a sorted index key, and then just sets it to this

      ! Input and output variable
      real(kind=8), dimension(:), intent(inout) :: list

      list = list(.sort.list)

      return
   end subroutine sort_double

   subroutine sort_string(list)

      ! Performs a heapsort into ascending order of the input list -- calls
      ! another routine to provide a sorted index key, and then just sets it to this

      ! Input and output variable
      character(len=*), dimension(:), intent(inout) :: list

      list = list(.sort.list)

      return
   end subroutine sort_string

   function sort_integer_return_key(list) result(key)

      ! Input/Output variables
      integer, dimension(:), intent(in) :: list
      integer, dimension(size(list))    :: key

      ! Local variables
      integer :: i, j, k, n, index, lists, keys

      n = size(list)

      ! initialize index into the original ordering
      key = (/(i, i=1, n)/)

      ! perform the heapsort of the input list
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            keys = key(k)
            lists = list(keys)
         else
            keys = key(index)
            lists = list(keys)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(key(j)) .lt. list(key(j + 1))) j = j + 1
            end if
            if (lists .lt. list(key(j))) then
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         key(i) = keys
      end do
      return
   end function sort_integer_return_key

   function sort_real_return_key(list) result(key)

      ! Input/Output variables
      real, dimension(:), intent(in) :: list
      integer, dimension(size(list)) :: key

      ! Local variables
      real    :: lists
      integer :: i, j, k, n, index, keys

      n = size(list)

      ! initialize index into the original ordering
      key = (/(i, i=1, n)/)

      ! perform the heapsort of the input list
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            keys = key(k)
            lists = list(keys)
         else
            keys = key(index)
            lists = list(keys)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(key(j)) .lt. list(key(j + 1))) j = j + 1
            end if
            if (lists .lt. list(key(j))) then
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         key(i) = keys
      end do
      return
   end function sort_real_return_key

   function sort_double_return_key(list) result(key)

      ! Input/Output variables
      real(kind=8), dimension(:), intent(in) :: list
      integer, dimension(size(list)) :: key

      ! Local variables
      real(kind=8)    :: lists
      integer :: i, j, k, n, index, keys

      n = size(list)

      ! initialize index into the original ordering
      key = (/(i, i=1, n)/)

      ! perform the heapsort of the input list
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            keys = key(k)
            lists = list(keys)
         else
            keys = key(index)
            lists = list(keys)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(key(j)) .lt. list(key(j + 1))) j = j + 1
            end if
            if (lists .lt. list(key(j))) then
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         key(i) = keys
      end do
      return
   end function sort_double_return_key

   function sort_strings_return_keys(list) result(key)

      ! Takes an input list of character strings and sorts it into
      ! alphabetical order using the Heapsort algorithm and returns
      ! a key (index) which can be used to produce an ordered array

      ! Input/Output variables
      character(len=*), dimension(:), intent(in) :: list
      integer, dimension(size(list))             :: key

      ! Local variables
      character(len=len(list)) :: lists
      integer :: i, j, k, n, index, keys

      n = size(list)

      ! initialize index into the original ordering
      key = (/(i, i=1, n)/)

      ! perform the heapsort of the input list
      k = n/2 + 1
      index = n
      do while (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            keys = key(k)
            lists = list(keys)
         else
            keys = key(index)
            lists = list(keys)
            key(index) = key(1)
            index = index - 1
            if (index .le. 1) then
               key(1) = keys
               return
            end if
         end if
         i = k
         j = k + k
         do while (j .le. index)
            if (j .lt. index) then
               if (list(key(j)) .lt. list(key(j + 1))) j = j + 1
            end if
            if (lists .lt. list(key(j))) then
               key(i) = key(j)
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         key(i) = keys
      end do
      return
   end function sort_strings_return_keys

   function unique_integers(list)

      ! Returns an list of the unique values in a sorted list. If
      ! the list isn't sorted this will return junk.

      ! Input variable
      integer, dimension(:), intent(in)             :: list

      ! Return value -- the dimension of the return value is
      ! determined by a call to another function, which determines
      ! the number of unique entries itself. This means we are
      ! doing the work twice, which is a pity. It lets us do neat
      ! tricks like this though, which is nice. If speed is an
      ! issue call .unique. operator and use the logical array mask
      ! directly
      integer, dimension(num_unique(list)) :: unique_integers

      ! Pack all the unique values from list into return array
      unique_integers = pack(list, (.unique. (list)))

   end function unique_integers

   pure function unique_integer_mask(list) result(mask)

      ! Returns an index of unique values in the sorted list. If
      ! the list isn't sorted this will return junk. The user will
      ! have to use a 'where (key /= 0)' or similar statement to get
      ! useful information from the key array.

      ! Input variables
      integer, dimension(:), intent(in) :: list

      ! Return value
      logical, dimension(size(list))    :: mask

      ! Test if each element in the array is the same as it's neighbour
      ! at an index value of one less
      mask = (list /= cshift(list, -1))

      ! Make sure the first value is always unique ...
      mask(1) = .TRUE.

   end function unique_integer_mask

   pure function num_unique_integers(list) result(n)

      ! Returns the number of unique values in a sorted list

      ! Input/Output variables
      integer, dimension(:), intent(in) :: list

      integer :: n

      ! Count the number of true values in the mask returned by .unique.
      n = count(.unique. (list))

   end function num_unique_integers

   function unique_reals(list)

      ! Returns an list of the unique values in a sorted list. If
      ! the list isn't sorted this will return junk.

      ! Input variable
      real, dimension(:), intent(in)             :: list

      ! Return value -- the dimension of the return value is
      ! determined by a call to another function, which determines
      ! the number of unique entries itself. This means we are
      ! doing the work twice, which is a pity. It lets us do neat
      ! tricks like this though, which is nice. If speed is an
      ! issue call .unique. operator and use the logical array mask
      ! directly
      real, dimension(num_unique(list)) :: unique_reals

      ! Pack all the unique values from list into return array
      unique_reals = pack(list, (.unique. (list)))

   end function unique_reals

   pure function unique_real_mask(list) result(mask)

      ! Returns an index of unique values in the sorted list. If
      ! the list isn't sorted this will return junk. The user will
      ! have to use a 'where (key /= 0)' or similar statement to get
      ! useful information from the key array.

      ! Input variables
      real, dimension(:), intent(in) :: list

      ! Return value
      logical, dimension(size(list))    :: mask

      ! Test if each element in the array is the same as it's neighbour
      ! at an index value of one less
      mask = (list /= cshift(list, -1))

      ! Make sure the first value is always unique ...
      mask(1) = .TRUE.

   end function unique_real_mask

   pure function num_unique_reals(list) result(n)

      ! Returns the number of unique values in a sorted list

      ! Input/Output variables
      real, dimension(:), intent(in) :: list

      integer :: n

      ! Count the number of true values in the mask returned by .unique.
      n = count(.unique. (list))

   end function num_unique_reals

   function unique_doubles(list)

      ! Returns an list of the unique values in a sorted list. If
      ! the list isn't sorted this will return junk.

      ! Input variable
      real(kind=8), dimension(:), intent(in)             :: list

      ! Return value -- the dimension of the return value is
      ! determined by a call to another function, which determines
      ! the number of unique entries itself. This means we are
      ! doing the work twice, which is a pity. It lets us do neat
      ! tricks like this though, which is nice. If speed is an
      ! issue call .unique. operator and use the logical array mask
      ! directly
      real(kind=8), dimension(num_unique(list)) :: unique_doubles

      ! Pack all the unique values from list into return array
      unique_doubles = pack(list, (.unique. (list)))

   end function unique_doubles

   pure function unique_double_mask(list) result(mask)

      ! Returns an index of unique values in the sorted list. If
      ! the list isn't sorted this will return junk. The user will
      ! have to use a 'where (key /= 0)' or similar statement to get
      ! useful information from the key array.

      ! Input variables
      real(kind=8), dimension(:), intent(in) :: list

      ! Return value
      logical, dimension(size(list))    :: mask

      ! Test if each element in the array is the same as it's neighbour
      ! at an index value of one less
      mask = (list /= cshift(list, -1))

      ! Make sure the first value is always unique ...
      mask(1) = .TRUE.

   end function unique_double_mask

   pure function num_unique_doubles(list) result(n)

      ! Returns the number of unique values in a sorted list

      ! Input/Output variables
      real(kind=8), dimension(:), intent(in) :: list

      integer :: n

      ! Count the number of true values in the mask returned by .unique.
      n = count(.unique. (list))

   end function num_unique_doubles

   function unique_strings(list)

      ! Returns an list of the unique values in a sorted list. If
      ! the list isn't sorted this will return junk.

      ! Input variable
      character(len=*), dimension(:), intent(in) :: list

      ! Return value -- the dimension of the return value is
      ! determined by a call to another function, which determines
      ! the number of unique entries itself. This means we are
      ! doing the work twice, which is a pity. It lets us do neat
      ! tricks like this though, which is nice. If speed is an
      ! issue call .unique. operator and use the logical array mask
      ! directly
      character(len=len(list)), dimension(num_unique(list)) :: unique_strings

      ! Pack all the unique values from list into return array
      unique_strings = pack(list, (.unique. (list)))

   end function unique_strings

   pure function unique_string_mask(list) result(mask)

      ! Returns an index of unique values in the sorted list. If
      ! the list isn't sorted this will return junk. The user will
      ! have to use a 'where (key /= 0)' or similar statement to get
      ! useful information from the key array.

      ! Input variables
      character(len=*), dimension(:), intent(in) :: list

      ! Return value
      logical, dimension(size(list))    :: mask

      ! Test if each element in the array is the same as it's neighbour
      ! at an index value of one less
      mask = (list /= cshift(list, -1))

      ! Make sure the first value is always unique ...
      mask(1) = .TRUE.

   end function unique_string_mask

   pure function num_unique_strings(list) result(n)

      ! Returns the number of unique values in a sorted list

      ! Input/Output variables
      character(len=*), dimension(:), intent(in) :: list

      integer :: n

      ! Count the number of true values in the mask returned by .unique.
      n = count(.unique. (list))

   end function num_unique_strings

   pure function reverse_integer(list) result(revlist)

      integer, dimension(:), intent(in)                 :: list
      integer, dimension(lbound(list, 1):ubound(list, 1)) :: revlist

      revlist = list(ubound(list, 1):lbound(list, 1):-1)

   end function reverse_integer

   pure function reverse_real(list) result(revlist)

      real, dimension(:), intent(in)                 :: list
      real, dimension(lbound(list, 1):ubound(list, 1)) :: revlist

      revlist = list(ubound(list, 1):lbound(list, 1):-1)

   end function reverse_real

   pure function reverse_double(list) result(revlist)

      real(kind=8), dimension(:), intent(in)                 :: list
      real(kind=8), dimension(lbound(list, 1):ubound(list, 1)) :: revlist

      revlist = list(ubound(list, 1):lbound(list, 1):-1)

   end function reverse_double

   pure function reverse_string(list) result(revlist)

      character(len=*), dimension(:), intent(in)                 :: list
      character(len=len(list)), dimension(lbound(list, 1):ubound(list, 1)) :: revlist

      revlist = list(ubound(list, 1):lbound(list, 1):-1)

   end function reverse_string

   RECURSIVE SUBROUTINE quick_sort(array, order)

      ! Quick sort routine from:
      ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
      ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
      ! Modified by Alan Miller to include an associated integer array which gives
      ! the positions of the elements in the original order.

      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)      :: array
      INTEGER, DIMENSION(:), INTENT(OUT)  :: order

      ! Local variable
      INTEGER :: i
      REAL, DIMENSION(size(array)) :: list

      DO i = 1, SIZE(list)
         order(i) = i
         list(i) = array(i)
      END DO

      CALL quick_sort_1(1, SIZE(list))

   CONTAINS

      RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

         INTEGER, INTENT(IN) :: left_end, right_end

         !     Local variables
         INTEGER             :: i, j, itemp
         REAL                :: reference, temp
         INTEGER, PARAMETER  :: max_simple_sort_size = 6

         IF (right_end < left_end + max_simple_sort_size) THEN
            ! Use interchange sort for small lists
            CALL interchange_sort(left_end, right_end)

         ELSE
            ! Use partition ("quick") sort
            reference = list((left_end + right_end)/2)
            i = left_end - 1; j = right_end + 1

            DO
               ! Scan list from left end until element >= reference is found
               DO
                  i = i + 1
                  IF (list(i) >= reference) EXIT
               END DO
               ! Scan list from right end until element <= reference is found
               DO
                  j = j - 1
                  IF (list(j) <= reference) EXIT
               END DO

               IF (i < j) THEN
                  ! Swap two out-of-order elements
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
               ELSE IF (i == j) THEN
                  i = i + 1
                  EXIT
               ELSE
                  EXIT
               END IF
            END DO

            IF (left_end < j) CALL quick_sort_1(left_end, j)
            IF (i < right_end) CALL quick_sort_1(i, right_end)
         END IF

      END SUBROUTINE quick_sort_1

      SUBROUTINE interchange_sort(left_end, right_end)

         INTEGER, INTENT(IN) :: left_end, right_end

         !     Local variables
         INTEGER             :: i, j, itemp
         REAL                :: temp

         DO i = left_end, right_end - 1
            DO j = i + 1, right_end
               IF (list(i) > list(j)) THEN
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
               END IF
            END DO
         END DO

      END SUBROUTINE interchange_sort

   END SUBROUTINE quick_sort

   function qsortd_integer(x) result(ind)

      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2002-12-18  Time: 11:55:47

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: x(:)
      INTEGER :: ind(size(x))
      INTEGER :: n

      !***************************************************************************

      !                                                         ROBERT RENKA
      !                                                 OAK RIDGE NATL. LAB.

      !   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
      ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
      ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
      ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
      ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
      ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
      ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
      ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
      ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
      ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
      ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
      ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
      ! UNSORTED PORTION.

      ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

      !                      X - VECTOR OF LENGTH N TO BE SORTED.

      !                    IND - VECTOR OF LENGTH >= N.

      ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

      ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
      !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
      !                          X IS DEFINED BY Y(I) = X(IND(I)).

      !*********************************************************************

      ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

      !*********************************************************************

      INTEGER   :: iu(21), il(21)
      INTEGER   :: m, i, j, k, l, ij, it, itt, indx
      REAL      :: r
      ! REAL (dp) :: t
      INTEGER :: t

      n = size(x)

      ! LOCAL PARAMETERS -

      ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
      !            INDICES OF PORTIONS OF THE ARRAY X
      ! M =      INDEX FOR IU AND IL
      ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
      ! K,L =    INDICES IN THE RANGE I,...,J
      ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
      ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
      ! INDX =   TEMPORARY INDEX FOR X
      ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
      ! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

      ! INITIALIZE IND, M, I, J, AND R

      DO i = 1, n
         ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

      ! TOP OF LOOP

20    IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
         r = r + .0390625
      ELSE
         r = r - .21875
      END IF

      ! INITIALIZE K

30    k = i

      ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j - i)
      it = ind(ij)
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) > t) THEN
         ind(ij) = indx
         ind(i) = it
         it = indx
         t = x(it)
      END IF

      ! INITIALIZE L

      l = j

      ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(j)
      IF (x(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = x(it)
      GO TO 50

      ! INTERCHANGE ELEMENTS K AND L

40    itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

      ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
      !   NOT LARGER THAN T

50    l = l - 1
      indx = ind(l)
      IF (x(indx) > t) GO TO 50

      ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60    k = k + 1
      indx = ind(k)
      IF (x(indx) < t) GO TO 60

      ! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

      ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
      !   ARRAY YET TO BE SORTED

      IF (l - i > j - k) THEN
         il(m) = i
         iu(m) = l
         i = k
         m = m + 1
         GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80

      ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70    m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

80    IF (j - i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

      ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90    i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i + 1)
      t = x(indx)
      it = indx
      indx = ind(i)
      IF (x(indx) <= t) GO TO 90
      k = i

100   ind(k + 1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < x(indx)) GO TO 100

      ind(k + 1) = it
      GO TO 90
   end function qsortd_integer

   function qsortd_real(x) result(ind)

      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2002-12-18  Time: 11:55:47

      IMPLICIT NONE
      ! INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
      INTEGER, PARAMETER  :: dp = rs_kind

      REAL(dp), INTENT(IN)  :: x(:)

      INTEGER :: ind(size(x))
      INTEGER :: n

      !***************************************************************************

      !                                                         ROBERT RENKA
      !                                                 OAK RIDGE NATL. LAB.

      !   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
      ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
      ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
      ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
      ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
      ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
      ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
      ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
      ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
      ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
      ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
      ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
      ! UNSORTED PORTION.

      ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

      !                      X - VECTOR OF LENGTH N TO BE SORTED.

      !                    IND - VECTOR OF LENGTH >= N.

      ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

      ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
      !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
      !                          X IS DEFINED BY Y(I) = X(IND(I)).

      !*********************************************************************

      ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

      !*********************************************************************

      INTEGER   :: iu(21), il(21)
      INTEGER   :: m, i, j, k, l, ij, it, itt, indx
      REAL      :: r
      REAL(dp) :: t

      n = size(x)

      ! LOCAL PARAMETERS -

      ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
      !            INDICES OF PORTIONS OF THE ARRAY X
      ! M =      INDEX FOR IU AND IL
      ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
      ! K,L =    INDICES IN THE RANGE I,...,J
      ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
      ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
      ! INDX =   TEMPORARY INDEX FOR X
      ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
      ! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

      ! INITIALIZE IND, M, I, J, AND R

      DO i = 1, n
         ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

      ! TOP OF LOOP

20    IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
         r = r + .0390625
      ELSE
         r = r - .21875
      END IF

      ! INITIALIZE K

30    k = i

      ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j - i)
      it = ind(ij)
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) > t) THEN
         ind(ij) = indx
         ind(i) = it
         it = indx
         t = x(it)
      END IF

      ! INITIALIZE L

      l = j

      ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(j)
      IF (x(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = x(it)
      GO TO 50

      ! INTERCHANGE ELEMENTS K AND L

40    itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

      ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
      !   NOT LARGER THAN T

50    l = l - 1
      indx = ind(l)
      IF (x(indx) > t) GO TO 50

      ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60    k = k + 1
      indx = ind(k)
      IF (x(indx) < t) GO TO 60

      ! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

      ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
      !   ARRAY YET TO BE SORTED

      IF (l - i > j - k) THEN
         il(m) = i
         iu(m) = l
         i = k
         m = m + 1
         GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80

      ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70    m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

80    IF (j - i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

      ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90    i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i + 1)
      t = x(indx)
      it = indx
      indx = ind(i)
      IF (x(indx) <= t) GO TO 90
      k = i

100   ind(k + 1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < x(indx)) GO TO 100

      ind(k + 1) = it
      GO TO 90
   END function qsortd_real

   function qsortd_double(x) result(ind)

      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2002-12-18  Time: 11:55:47

      IMPLICIT NONE
      ! INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
      INTEGER, PARAMETER  :: dp = rd_kind

      REAL(dp), INTENT(IN)  :: x(:)

      INTEGER :: ind(size(x))
      INTEGER :: n

      !***************************************************************************

      !                                                         ROBERT RENKA
      !                                                 OAK RIDGE NATL. LAB.

      !   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
      ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
      ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
      ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
      ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
      ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
      ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
      ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
      ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
      ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
      ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
      ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
      ! UNSORTED PORTION.

      ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

      !                      X - VECTOR OF LENGTH N TO BE SORTED.

      !                    IND - VECTOR OF LENGTH >= N.

      ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

      ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
      !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
      !                          X IS DEFINED BY Y(I) = X(IND(I)).

      !*********************************************************************

      ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

      !*********************************************************************

      INTEGER   :: iu(21), il(21)
      INTEGER   :: m, i, j, k, l, ij, it, itt, indx
      REAL      :: r
      REAL(dp) :: t

      n = size(x)

      ! LOCAL PARAMETERS -

      ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
      !            INDICES OF PORTIONS OF THE ARRAY X
      ! M =      INDEX FOR IU AND IL
      ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
      ! K,L =    INDICES IN THE RANGE I,...,J
      ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
      ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
      ! INDX =   TEMPORARY INDEX FOR X
      ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
      ! T =      CENTRAL ELEMENT OF X

      IF (n <= 0) RETURN

      ! INITIALIZE IND, M, I, J, AND R

      DO i = 1, n
         ind(i) = i
      END DO
      m = 1
      i = 1
      j = n
      r = .375

      ! TOP OF LOOP

20    IF (i >= j) GO TO 70
      IF (r <= .5898437) THEN
         r = r + .0390625
      ELSE
         r = r - .21875
      END IF

      ! INITIALIZE K

30    k = i

      ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

      ij = i + r*(j - i)
      it = ind(ij)
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) > t) THEN
         ind(ij) = indx
         ind(i) = it
         it = indx
         t = x(it)
      END IF

      ! INITIALIZE L

      l = j

      ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(j)
      IF (x(indx) >= t) GO TO 50
      ind(ij) = indx
      ind(j) = it
      it = indx
      t = x(it)

      ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
      !   INTERCHANGE IT WITH T

      indx = ind(i)
      IF (x(indx) <= t) GO TO 50
      ind(ij) = indx
      ind(i) = it
      it = indx
      t = x(it)
      GO TO 50

      ! INTERCHANGE ELEMENTS K AND L

40    itt = ind(l)
      ind(l) = ind(k)
      ind(k) = itt

      ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
      !   NOT LARGER THAN T

50    l = l - 1
      indx = ind(l)
      IF (x(indx) > t) GO TO 50

      ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60    k = k + 1
      indx = ind(k)
      IF (x(indx) < t) GO TO 60

      ! IF K <= L, INTERCHANGE ELEMENTS K AND L

      IF (k <= l) GO TO 40

      ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
      !   ARRAY YET TO BE SORTED

      IF (l - i > j - k) THEN
         il(m) = i
         iu(m) = l
         i = k
         m = m + 1
         GO TO 80
      END IF

      il(m) = k
      iu(m) = j
      j = l
      m = m + 1
      GO TO 80

      ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70    m = m - 1
      IF (m == 0) RETURN
      i = il(m)
      j = iu(m)

80    IF (j - i >= 11) GO TO 30
      IF (i == 1) GO TO 20
      i = i - 1

      ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90    i = i + 1
      IF (i == j) GO TO 70
      indx = ind(i + 1)
      t = x(indx)
      it = indx
      indx = ind(i)
      IF (x(indx) <= t) GO TO 90
      k = i

100   ind(k + 1) = ind(k)
      k = k - 1
      indx = ind(k)
      IF (t < x(indx)) GO TO 100

      ind(k + 1) = it
      GO TO 90
   END function qsortd_double

end module sort_functions

module mod_sort

   use sort_functions

end module mod_sort
