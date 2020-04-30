module SWIO

  !-----------------------------------------------------------------------------
  ! SWIO Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS           => SetServices,          &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance,        &
    model_label_Finalize       => label_Finalize,       &
    model_label_CheckImport    => label_CheckImport

  use COMIO

  use swio_data
  use swio_methods

  implicit none

  private

  public :: SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! begin
    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! switch to IPDv03
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry points for methods that require specific implementation
    ! - advertise import fields and set TransferOfferGeomObject attribute.
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! - realize "connected" import fields with TransferActionGeomObject
    !   equal to "provide".
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeP3, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return

    ! - realize "connected" import fields with TransferActionGeomObject
    !   equal to "accept" on the transferred Grid/Mesh/LocStream objects
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeP5, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! - set component's initialization status
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=DataInitialize, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    ! - advance method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! - finalize method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! - do not check time stamp of imported fields
    call ESMF_MethodRemove(gcomp, model_label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
      specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name, value
    type(ESMF_Config)          :: config

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP0"
    
    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get name of config file
    call ESMF_AttributeGet(gcomp, name="ConfigFile", value=value, &
      defaultValue="swio.conf", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (btest(verbosity,8)) then
      call ESMF_LogWrite(trim(name)//": ConfigFile = "//trim(value), &
        ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! load component's configuration
    config = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ConfigLoadFile(config, value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! store Config object into gridded component
    call ESMF_GridCompSet(gcomp, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    logical                    :: listEnd
    integer                    :: item
    integer                    :: stat
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name
    character(len=ESMF_MAXSTR) :: fieldStandardName
    character(len=ESMF_MAXSTR) :: transferOfferGeomObject, sharePolicyField
    character(len=ESMF_MAXSTR) :: msgString
    type(ESMF_Config)          :: config
    type(SWIO_InternalState_T) :: is
    type(SWIO_Data_T), pointer :: this

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP1"

    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get component's configuration
    call ESMF_GridCompGet(gcomp, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! allocate memory for the internal state and store it into component
    allocate(is % wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % wrap
    ! initialize internal state
    this % fieldCount = 0
    this % gridType   = ""
    this % filePrefix = ""
    this % fileSuffix = ""
    nullify(this % io)

    ! get output grid selection
    call ESMF_ConfigGetAttribute(config, this % gridType, &
      label="output_grid_type:", default="none", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! check if output grid type is supported
    select case (trim(this % gridType))
      case ("none")
      case ("latlon")
      case default
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg=" -  grid type: "//trim(this % gridType), &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return
    end select
    if (btest(verbosity,8)) then
      call ESMF_LogWrite(trim(name)//": gridType = "//trim(this % gridType), &
        ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! set import fields attributes based on grid type
    ! -- advertise field
    if (trim(this % gridType) == "none") then
      transferOfferGeomObject="cannot provide"
      sharePolicyField="share"
    else
      transferOfferGeomObject="will provide"
      sharePolicyField="not share"
    end if

    ! get import fields from Config
    ! - locate field table
    call ESMF_ConfigFindLabel(config, "import_fields::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! - read field standard name from 1st column and advertise field
    item = 0
    do
      ! get next row
      call ESMF_ConfigNextLine(config, tableEnd=listEnd, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      if (listEnd) exit
      ! get input field name
      call ESMF_ConfigGetAttribute(config, fieldStandardName, rc=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ! advertise field
      call NUOPC_Advertise(importState, StandardName=fieldStandardName, &
        SharePolicyField=sharePolicyField, &
        TransferOfferGeomObject=transferOfferGeomObject, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      item = item + 1
      if (btest(verbosity,8)) then
        write(msgString,'(a,": import[",i0,"]: ",a)') trim(name), &
          item, trim(fieldStandardName)
        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(trim(name)//": - GeomObject: "&
          //trim(transferOfferGeomObject), ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(trim(name)//": - Field: "&
          //trim(sharePolicyField), ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end do
    if (btest(verbosity,8)) then
      if (item == 0) then
        call ESMF_LogWrite(trim(name)//": import: None", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP1

  subroutine InitializeP3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    logical                         :: createGrid, isPresent
    integer                         :: item, itemCount
    integer                         :: stat
    integer                         :: verbosity
    integer                         :: fieldCount, nameCount
    integer                         :: ungriddedDimLength
    character(ESMF_MAXSTR)          :: name
    character(ESMF_MAXSTR)          :: msgString
    character(ESMF_MAXSTR)          :: transferAction
    character(ESMF_MAXSTR), pointer :: connectedList(:)
    character(ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Field)                :: field
    type(ESMF_Grid)                 :: grid
    type(SWIO_InternalState_T)      :: is
    type(SWIO_Data_T), pointer      :: this

    ! local parameters
    character(len=*), parameter :: rName = "InitializeP3"

    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

     ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % wrap

    ! determine how many fields are connected
    nullify(standardNameList)
    nullify(connectedList)
    call NUOPC_GetStateMemberLists(importState, &
      StandardNameList=standardNameList, ConnectedList=connectedList, &
      nestedFlag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    nameCount = 0
    if (associated(standardNameList)) nameCount = size(standardNameList)
    fieldCount = 0
    if (nameCount > 0) then
      if (associated(connectedList)) then
        do item = 1, nameCount
          if (trim(connectedList(item)) == "true") fieldCount = fieldCount + 1
        end do
      end if
      if (fieldCount == 0) nameCount = 0
    end if

    fieldCount = 0
    ungriddedDimLength = 0
    createGrid = .true.
    do item = 1, nameCount
      if (trim(connectedList(item)) == "true") then
        call ESMF_StateGet(importState, &
          itemName=trim(standardNameList(item)), field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
          value=transferAction, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (btest(verbosity,8)) then
          call ESMF_LogWrite(trim(name)//": "//rName &
            //": TransferActionGeomObject = "//trim(transferAction), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
        end if
        if (trim(transferAction) == "provide") then
          select case (this % gridType)
            case ("none")
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              exit
            case ("latlon")
              if (createGrid) then
                grid = SWIO_GridCreateLatLon(gcomp, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__,  &
                  file=__FILE__)) &
                  return  ! bail out
                call ESMF_AttributeGet(grid, name="UngriddedDimLength", &
                  convention="NUOPC", purpose="Instance", isPresent=isPresent, &
                  rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__,  &
                  file=__FILE__)) &
                  return  ! bail out
                if (isPresent) then
                  call ESMF_AttributeGet(grid, name="UngriddedDimLength", &
                    value=ungriddedDimLength, convention="NUOPC", &
                    purpose="Instance", rc=rc)
                  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                    line=__LINE__,  &
                    file=__FILE__)) &
                    return  ! bail out
                end if
                createGrid = .false.
              end if
              if (ungriddedDimLength > 0) then
                field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
                  ungriddedLBound=(/ 1 /), &
                  ungriddedUBound=(/ ungriddedDimLength /), &
                  name=trim(standardNameList(item)), rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__,  &
                  file=__FILE__)) &
                  return  ! bail out
              else
                field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
                  name=trim(standardNameList(item)), rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__,  &
                  file=__FILE__)) &
                  return  ! bail out
              end if
            case default
              exit
          end select
          call NUOPC_Realize(importState, field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
          fieldCount = fieldCount + 1
          if (btest(verbosity,8)) then
            write(msgString,'(a,": import[",i0,"]:",2(1x,a))') trim(name), &
              item, trim(standardNameList(item)), "(realized on provided GeomObject)"
            call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
          end if
        end if
      else
        call ESMF_StateRemove(importState, standardNameList(item), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end do

    this % fieldCount = fieldCount

    if (associated(standardNameList)) then
      deallocate(standardNameList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if
    if (associated(connectedList)) then
      deallocate(connectedList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP3
    
  subroutine InitializeP5(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                             :: item
    integer                             :: stat
    integer                             :: verbosity
    integer                             :: nameCount
    integer                             :: fieldCount
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: transferAction
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Field)                    :: field
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! local parammeters
    character(len=*), parameter :: rName = "InitializeP5"

    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % wrap

    fieldCount = 0
    if (trim(this % gridType) == "none") then
      ! realize fields with TransferActionGeomObject "accept"
      nullify(standardNameList)
      call NUOPC_GetStateMemberLists(importState, &
        StandardNameList=standardNameList, nestedFlag=.true., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      nameCount = 0
      if (associated(standardNameList)) nameCount = size(standardNameList)

      do item = 1, nameCount
        call ESMF_StateGet(importState, &
          itemName=trim(standardNameList(item)), field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
          value=transferAction, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (btest(verbosity,8)) then
          call ESMF_LogWrite(trim(name)//": "//rName &
            //": TransferActionGeomObject = "//trim(transferAction), ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
        end if
        select case (trim(transferAction))
!         case ("accept")
          case ("accept","provide")
            call NUOPC_Realize(importState, fieldName=standardNameList(item), rc=rc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
            if (btest(verbosity,8)) then
              write(msgString,'(a,": import[",i0,"]:",2(1x,a))') trim(name), &
                item, trim(standardNameList(item)), "(realized on accepted GeomObject)"
              call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__)) &
                return  ! bail out
            end if
            fieldCount = fieldCount + 1
          case ("complete")
            if (btest(verbosity,8)) then
              write(msgString,'(a,": import[",i0,"]:",2(1x,a))') trim(name), &
                item, trim(standardNameList(item)), "(completed)"
              call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__)) &
                return  ! bail out
            end if
            fieldCount = fieldCount + 1
          case default
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
              msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            exit
        end select
      end do
      if (associated(standardNameList)) then
        deallocate(standardNameList, stat=stat)
        if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
          msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
    end if

    if (fieldCount > 0) then
      if (this % fieldCount > 0) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_INCONS, &
          msg="Fields that can both accept and receive GeomObjects found", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      else
        this % fieldCount = fieldCount
      end if
    end if

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP5
  
  !-----------------------------------------------------------------------------

  subroutine DataInitialize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                    :: comm
    integer                    :: iofmt
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: ioFormat
    character(len=ESMF_MAXSTR) :: name, svalue
    type(ESMF_Config)          :: config
    type(ESMF_VM)              :: vm
    type(SWIO_InternalState_T) :: is
    type(SWIO_Data_T), pointer :: this
    
    ! local parammeters
    character(len=*), parameter :: rName = "DataInitialize"
    
    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    this => is % wrap

    ! initialize I/O layer
    if (this % fieldCount > 0) then
      ! get component's configuration
      call ESMF_GridCompGet(gcomp, config=config, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      ! get I/O format
      call ESMF_ConfigGetAttribute(config, svalue, &
        label="output_format:", default="none", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ioFormat = ESMF_UtilStringLowerCase(svalue, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      ! convert text to COMIO I/O selector
      select case (trim(ioFormat))
        case ("hdf5")
          iofmt = COMIO_FMT_HDF5
          this % fileSuffix = "hd5"
        case ("pnetcdf", "parallel-netcdf")
          iofmt = COMIO_FMT_PNETCDF
          this % fileSuffix = "nc"
        case default
          call ESMF_LogSetError(rcToCheck=ESMF_RC_VAL_OUTOFRANGE, &
            msg="Invalid I/O format", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
      end select

      ! get file info
      call ESMF_ConfigGetAttribute(config, this % filePrefix, &
        label="output_file_prefix:", default="data", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      ! get MPI communicator from component's VM
      call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
  
      ! initialize COMIO
      this % io => COMIO_T(fmt=iofmt, comm=comm)
      if (this % io % err % check(msg="Failure initializing I/O", &
        line=__LINE__,  &
        file=__FILE__)) then
        call ESMF_LogSetError(ESMF_RC_OBJ_INIT, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)
        return  ! bail out
      end if
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//rName//": I/O initialized (COMIO)"&
          //" - Writing to: "//trim(this % filePrefix)//"_YYYYMMDD_hhmmss."&
          //trim(this % fileSuffix)//" "//trim(ioFormat)//" files", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    else
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//rName//": I/O not initialized"&
          //" - No fields present", ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
        
  end subroutine DataInitialize

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    logical                             :: isTimeValid
    integer                             :: i, item, stat
    integer                             :: dimCount
    integer                             :: nameCount
    integer                             :: rank
    integer                             :: verbosity
    character(len=15)                   :: timeStamp
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: fileName
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Array)                    :: array
    type(ESMF_Field)                    :: field
    type(ESMF_GeomType_Flag)            :: geomtype
    type(ESMF_Grid)                     :: grid
    type(ESMF_Mesh)                     :: mesh
    type(ESMF_MeshLoc)                  :: meshloc
    type(ESMF_State)                    :: importState
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! local parameters
    character(len=*), parameter :: rName = "Run"

    ! begin
    rc = ESMF_SUCCESS

    ! write output
    call SWIO_Output(gcomp, phaseName=rName, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelAdvance
 
  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    logical                    :: configIsPresent
    integer                    :: stat
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name
    type(ESMF_Config)          :: config
    type(SWIO_InternalState_T) :: is

    ! local parameters
    character(len=*), parameter :: rName = "Finalize"
    ! begin
    rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! write output at the final time step
    call SWIO_Output(gcomp, phaseName=rName, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! check if Config is present
    call ESMF_GridCompGet(gcomp, configIsPresent=configIsPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (configIsPresent) then
      ! get component's Config object
      call ESMF_GridCompGet(gcomp, config=config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ! destroy config
      call ESMF_ConfigDestroy(config, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//rName//": Config successfully destroyed",&
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

    ! get internal state
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (associated(is % wrap)) then
      if (associated(is % wrap % io)) then
        ! shutdown I/O layer
        call is % wrap % io % shutdown()
        if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        ! finally, deallocate internal state
        deallocate(is % wrap % io, stat=stat)
        if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
          msg="Failed to free memory used by I/O layer.", &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        end if
        if (btest(verbosity,8)) then
          call ESMF_LogWrite(trim(name)//": "//rName//": I/O successfully shutdown",&
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
        end if
      ! finally, deallocate internal state
      deallocate(is % wrap, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Failed to free internal state memory.", &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      nullify(is % wrap)
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//rName//": Internal state memory successfully released",&
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

  end subroutine Finalize

end module SWIO
