module swio_data

  use comio
  use ESMF

  implicit none

  type SWIO_Data_T
    integer                 :: logLevel
    integer                 :: fieldCount
    character(ESMF_MAXSTR)  :: gridType
    character(ESMF_MAXSTR)  :: filePrefix
    character(ESMF_MAXSTR)  :: fileSuffix
    class(COMIO_T), pointer :: io
  end type

  type SWIO_InternalState_T
    type(SWIO_Data_T), pointer :: wrap
  end type

end module swio_data
