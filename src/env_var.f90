!! SPDX-License-Identifier: BSD-3-Clause

! Module for reading environment variables

module m_env_var

   use decomp_2d_mpi

   implicit none

   ! Default is private
   private

   ! Default value if no default value is provided
   integer, parameter :: env_var_default_int = 0

   public :: get_env_var

   interface get_env_var
      module procedure get_env_var_char
      module procedure get_env_var_int
   end interface get_env_var

contains

   ! Extract the raw environment variable
   function get_env_var_char(name) result(output)

      implicit none

      ! Argument and result
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: output

      ! Local variables
      integer :: len, ierror

      ! Get the size of the environment variable
      call get_environment_variable(name, length=len, &
                                    status=ierror, trim_name=.true.)
      if (ierror == 1) then
         output = ''
         call decomp_2d_warning(ierror, "Environment variable " // &
                                        trim(name) //" is not defined.")
         return
      else if (ierror /= 0) then
         output = ''
         call decomp_2d_warning(ierror, "No support for environment variable")
         return
      end if

      ! Allocate the output variable
      allocate(character(len=max(1, len)) :: output)

      ! Get the environment variable
      call get_environment_variable(name, output, &
                                    status=ierror, trim_name=.true.)

      if (ierror /= 0) then
         output = ''
         call decomp_2d_warning(__FILE__, __LINE__, ierror, &
                                "Error when reading " // name)
         return
      end if

   end function get_env_var_char

   ! Extract the raw environment variable and convert to int
   function get_env_var_int(name, default) result(output)

      implicit none

      ! Arguments and result
      character(len=*), intent(in) :: name
      integer, intent(in) :: default
      integer :: output

      ! Local variables
      integer :: ierror
      character(len=:), allocatable :: raw, fmt

      ! Get the raw output
      raw = get_env_var_char(name)

      ! Return if the environment variable was not available
      if (raw == '') then
         output = default
         deallocate(raw)
         return
      end if

      ! Format for the conversion
      allocate(character(len = 3 + max(1,len(raw))) :: fmt)
      fmt(1:2) = '(i'
      write(fmt(3:len(fmt)-1), *) max(1,len(raw))
      fmt(len(fmt):len(fmt)) = ')'

      ! Integer convertion
      read(raw, fmt, iostat=ierror) output

      ! Free memory
      deallocate(fmt)
      deallocate(raw)

   end function get_env_var_int

end module m_env_var
