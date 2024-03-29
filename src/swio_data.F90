module swio_data

  use comio
  use ESMF

  implicit none

  type SWIO_Task_T
    type(ESMF_Field),   pointer :: fieldOut(:)
    type(ESMF_Field),   pointer :: fieldInp(:)
    real(ESMF_KIND_R8), pointer :: paramInp(:)
    character(ESMF_MAXSTR)      :: operation
    real(ESMF_KIND_R8)          :: scaleFactor
  end type

  type SWIO_Mask_T
    type(ESMF_Field)   :: field
    real(ESMF_KIND_R8) :: value
    real(ESMF_KIND_R8) :: fill
  end type

  type SWIO_Pair_T
    character(ESMF_MAXSTR) :: key
    character(ESMF_MAXSTR) :: value
  end type

  type SWIO_Data_T
    integer                 :: fieldCount
    integer                 :: outputCount
    logical                 :: geoReference
    character(ESMF_MAXSTR)  :: gridType
    character(ESMF_MAXSTR)  :: filePrefix
    character(ESMF_MAXSTR)  :: fileSuffix
    character(2)            :: cmode
    type(SWIO_Mask_T), pointer     :: mask
    type(SWIO_Pair_T), pointer     :: meta(:)
    type(SWIO_Pair_T), pointer     :: output(:)
    type(SWIO_Task_T), pointer     :: task(:)
    class(COMIO_T),    allocatable :: io
  end type

  type SWIO_InternalState_T
    type(SWIO_Data_T), pointer :: wrap
  end type

end module swio_data
