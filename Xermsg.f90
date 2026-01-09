!*****************************************************************************!
!
Module ERRORS
!
!*****************************************************************************!
!
!--- FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!
use DImach
implicit none

contains
!=============================================================================!
!
!DECK XERMSG
SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!***BEGIN PROLOGUE  XERMSG
!***PURPOSE  Process error messages for SLATEC and other libraries.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
      INTEGER :: LKNTRL, I, MAXMES, KDUMMY, KOUNT, LERR, LLEVEL, MKNTRL, LTEMP, LEVEL, NERR
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR. &
     &   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
     &      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
     &      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
!
!       RECORD THE MESSAGE.
!
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!       IS NOT ZERO.
!
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!
      IF (LKNTRL .GT. 0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR. &
     &       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
!
!       NOW SEND OUT THE MESSAGE.
!
      CALL XERPRN (' *  ', -1, MESSG, 72)
!
!       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
!
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
!
!       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
!
!       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
END subroutine XERMSG

!=============================================================================!
!
!DECK FDUMP
SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
END SUBROUTINE FDUMP

!=============================================================================!
!
!*DECK J4SAVE
FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
!***BEGIN PROLOGUE  J4SAVE
!***SUBSIDIARY
!***PURPOSE  Save or recall global variables needed by error
!            handling routines.
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER :: J4SAVE, IWHICH, IVALUE
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
END FUNCTION J4SAVE


!=============================================================================!
!
!*DECK XERCNT
SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!***BEGIN PROLOGUE  XERCNT
!***SUBSIDIARY
!***PURPOSE  Allow user control over handling of errors.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCNT-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCNT
      INTEGER :: NERR, LEVEL, KONTRL
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
!***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
END SUBROUTINE XERCNT


!=============================================================================!
!
!*DECK XERHLT
SUBROUTINE XERHLT (MESSG)
!***BEGIN PROLOGUE  XERHLT
!***SUBSIDIARY
!***PURPOSE  Abort program execution and print error message.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERHLT-A)
!***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG is as in XERMSG.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
END SUBROUTINE XERHLT

!=============================================================================!
!
!*DECK XERPRN
SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!***BEGIN PROLOGUE  XERPRN
!***SUBSIDIARY
!***PURPOSE  Print error messages processed by XERMSG.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERPRN-A)
!***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
! This routine sends one or more lines to each of the (up to five)
! logical units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP, I, IDELTA, LWRAP, NEXTC, LPIECE, LENMSG, LPREF, N
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
!***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
!
!       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
!
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
!
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
!
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
!
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
!
      LWRAP = MAX(16, MIN(132, NWRAP))
!
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
!
!       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
!
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
!
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
!
!       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
!
!       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
!
!       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
!
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
!
!       THERE WAS NO NEW LINE SENTINEL FOUND.
!
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
!
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
!
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
!
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
!
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
!
!       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
!
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
!
!       PRINT
!
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
!
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
END SUBROUTINE XERPRN

!=============================================================================!
!
!*DECK XERSVE
SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!***BEGIN PROLOGUE  XERSVE
!***SUBSIDIARY
!***PURPOSE  Record that an error has occurred.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
      INTEGER, PARAMETER :: LENTAB=10
      INTEGER LUN(5), KOUNTX
      INTEGER I, NMSG, NERTAB, LEVTAB, KOUNT, IUNIT,KUNIT,NUNIT, ICOUNT, LEVEL, NERR, KFLAG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
      IF (KFLAG.LE.0) THEN
!
!        Dump the table.
!
         IF (NMSG.EQ.0) RETURN
!
!        Print to each unit.
!
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
            WRITE (IUNIT,9000)
!
!           Print body of table.
!
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
     &            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND. &
     &         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND. &
     &         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
!
         IF (NMSG.LT.LENTAB) THEN
!
!           Empty slot found for new message.
!
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
!
!           Table is full.
!
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
     &   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',&
     &   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
END SUBROUTINE XERSVE

!=============================================================================!
!
!*DECK XGETUA
SUBROUTINE XGETUA (IUNITA, N)
!***BEGIN PROLOGUE  XGETUA
!***PURPOSE  Return unit number(s) to which error messages are being
!            sent.
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETUA
      INTEGER :: I, N, INDEX
      INTEGER, DIMENSION(5) :: IUNITA
!***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
END SUBROUTINE XGETUA

end module ERRORS

!=============================================================================!
!
!
