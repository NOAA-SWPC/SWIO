module swio_data

  use comio
  use ESMF

  implicit none

  type SWIO_Task_T
    type(ESMF_Field), pointer :: fieldOut(:)
    type(ESMF_Field), pointer :: fieldInp(:)
    character(ESMF_MAXSTR)    :: operation
    real(ESMF_KIND_R8)        :: scaleFactor
  end type

  type SWIO_Data_T
    integer                 :: logLevel
    integer                 :: fieldCount
    logical                 :: geoReference
    character(ESMF_MAXSTR)  :: gridType
    character(ESMF_MAXSTR)  :: filePrefix
    character(ESMF_MAXSTR)  :: fileSuffix
    type(SWIO_Task_T), pointer :: task(:)
    class(COMIO_T),    pointer :: io
  end type

  type SWIO_InternalState_T
    type(SWIO_Data_T), pointer :: wrap
  end type

end module swio_data
