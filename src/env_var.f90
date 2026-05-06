!! SPDX-License-Identifier: BSD-3-Clause

! Module for reading environment variables

submodule(decomp_2d) m_env_var

   use decomp_2d_mpi

   implicit none

contains

   ! Extract the raw environment variable
   module function get_env_var_char(name) result(output)

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

         ! The requested environment variable was not found
         ! Set default output
         output = ''
         ! Print a warning if the debug level is high enough
         if (decomp_debug >= D2D_DEBUG_LEVEL_TRACE) &
            call decomp_2d_warning(ierror, "Environment variable "// &
                                   trim(name)//" is not defined.")
         return

      else if (ierror /= 0) then

         ! This is a low probability event
         output = ''
         call decomp_2d_warning(ierror, "No support for environment variable")
         return

      end if

      ! Allocate the output variable
      allocate (character(len=max(1, len)) :: output)

      ! Get the environment variable
      call get_environment_variable(name, output, &
                                    status=ierror, trim_name=.true.)

      if (ierror /= 0) then

         ! This is a low probability event
         output = ''
         call decomp_2d_warning(__FILE__, __LINE__, ierror, &
                                "Error when reading "//trim(name))
         return

      end if

   end function get_env_var_char

   ! Extract the raw environment variable and convert to int
   module function get_env_var_int(name, default) result(output)

      implicit none

      ! Arguments and result
      character(len=*), intent(in) :: name
      integer, intent(in) :: default
      integer :: output

      ! Local variables
      integer :: ierror, io_unit
      character(len=:), allocatable :: raw, fmt

      ! Get the raw output
      raw = get_env_var_char(name)

      ! Return if the environment variable was not available
      if (raw == '') then

         ! Update the log if the debug level is high enough
         if (decomp_debug >= D2D_DEBUG_LEVEL_TRACE) then
            if (d2d_log_is_active()) then
               io_unit = d2d_log_get_unit()
               write (io_unit, *) trim(name)//" set to default ", output
               call d2d_log_close_unit(io_unit)
            end if
         end if

         ! Set the output to the default value
         output = default

         ! Free memory and return
         deallocate (raw)
         return

      end if

      ! Format for the conversion
      allocate (character(len=3 + max(1, len(raw))) :: fmt)
      write (fmt, '(*(g0))') '(i', max(1, len(raw)), ')'

      ! Character => Integer
      read (raw, fmt, iostat=ierror) output

      ! Update the log if the debug level is high enough
      if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) then
         if (d2d_log_is_active()) then
            io_unit = d2d_log_get_unit()
            write (io_unit, *) trim(name)//" set to ", output
            call d2d_log_close_unit(io_unit)
         end if
      end if

      ! Free memory
      deallocate (fmt)
      deallocate (raw)

   end function get_env_var_int

end submodule m_env_var
