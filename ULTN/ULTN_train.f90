! PROGRAM TO TRAIN A Matrix Product States (MPS) OBJECT
!	ON A FRACTION OF MNIST DATASET (1000 SAMPLES) TO TRY 
!	TO LEARN PIXEL DISTRIBUTION FEATURES (UNSUPERVISED LEARNING TASK)

PROGRAM ULTN_TRAIN

! Importing useful modules
USE UTILITY

IMPLICIT NONE

!-----VARIABLES-----
! dataset_load		: the set of images (70k * (28 * 28) array, with elements {0,1})
!						It is used to easily load the file, then train elements
!						are selected and used in "dataset" variable
! dataset			: the selected elements to train on
! dataset_mask_r	: the indexes mask of random numbers
! dataset_mask_i	: the indexes mask of integer numbers to select elements from dataset_load
! params			: parameters for the analysis/training (see later or mod_utility file)
! b_v				: array representing the bonds values (dimension = 785, because it contains
!						the dimension of the left and right bonds of 784 middle tensors)
! max_b_v			: maximum bonds allowed (otherwise the tensors can not be reconstructed
!							after being merged)
! matrices			: list of tensors of the network
! cumuls			: array of cumulants
! ii, jj, kk		: indexes to be used in loops
! chosen_seed		: array of number to use as seed to have reproducible executions
! num_file			: string useful to save several files

!-----DATASET VARIABLES-----
INTEGER(KIND=1), DIMENSION(70000,784) :: dataset_load
INTEGER(KIND=1), DIMENSION(1000,784) :: dataset
REAL*8, DIMENSION(1000) :: dataset_mask_r
INTEGER, DIMENSION(1000) :: dataset_mask_i
!-----MPS & TRAIN VARIABLES-----
TYPE(PARAMETERS) :: params
INTEGER, DIMENSION(785) :: b_v, max_b_v
TYPE(VEC_TENS) :: matrices
TYPE(CUMULANTS) :: cumuls
!-----OTHER VARIABLES-----
INTEGER :: ii, jj, kk
INTEGER, DIMENSION(33) :: chosen_seed	! N.B.: usually the seed is an array of size 33
CHARACTER(LEN=3) :: num_file

! STARTING COMPUTATIONS

! Initializing seed using random but fixed values (priviously sampled)
chosen_seed = (/70290, 38879, 10013, 13932, 41407, 69718, 84713, 49520, 32977, 46194, 81813, 39405, &
& 85864, 19900, 73449, 56005, 19063, 74073, 25801, 63107, 71317, 22322, 90386, 40631, 99486, 65006, &
& 89833, 49354, 69803, 70049, 54845, 52086, 76196/)
CALL RANDOM_SEED(put=chosen_seed)

! Reading dataset from file
OPEN(33, FILE="data.bin",FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
READ(33) dataset_load
WRITE(*,*) SHAPE(dataset)
CLOSE(33)

! Initializing 1000 random numbers in [0,1)
CALL RANDOM_NUMBER(dataset_mask_r)

! Scaling [0,1) to [0, 70k)
!	Truncating to integers
!	Adding 1 to all indexes
dataset_mask_i = 1+FLOOR(dataset_mask_r*70000)

! Selecting 1000 elements using mask
dataset = dataset_load(dataset_mask_i,:)

! Verbose
WRITE(*,*) "Training set is loaded."

! Plotting some images from the dataset to see if they are loaded correctly
!WRITE(6,*) "TEST DATASET"
!DO kk=1,10
!	WRITE(6,*) "Img", kk
!	DO ii=1,28
!		DO jj = 1,27
!			WRITE(6,"(I2.1)", ADVANCE="No") dataset(kk,(ii-1)*28+jj)
!		END DO
!		WRITE(6,"(I2.1)", ADVANCE="Yes") dataset(kk,ii*28)
!	END DO
!END DO


! INITIALIZING PARAMS VARIABLES
! space_size		: number of input size (28x28)
! descent_steps		: number of steps for GD process (on a single merged tensor)
! init_bond_dim		: starting dimension of each bond
! maxi_bond			: maximum bond dimension allowed
! mini_bond			: minimum bond dimension allowed
! cut_off			: threshold for SVD procedure
! learning_rate		: gradient descent learning rate

params%space_size = 784
params%descent_steps = 10
params%num_of_iter = 50
params%init_bond_dim = 2
params%mini_bond = 2
params%maxi_bond = 20
params%cut_off = 1d-7
params%learning_rate = 0.01


! ----- INITIALIZING b_v -----
b_v = params%init_bond_dim
b_v(1) = 1
b_v(785) = 1


! ----- INITIALIZING max_b_v (useful for left canonization checks) -----
! The maximum dimensions allowed for the n-th position are: min(input_dimension^(n-1), maxi_bond)
! Initializing max_b_v (untill the values become greater than maxi_bond)
ii = 1
! In the next line, instead of 2 the dimensio of the input 
!	should be inserted (our MPS tensors have dimensions [left, 2, right]
DO WHILE (2**(ii-1).LE.params%maxi_bond)
	max_b_v(ii) = 2**(ii-1)
	max_b_v(785+1-ii) = 2**(ii-1)
	ii = ii+1
END DO

! Setting the others to maxi_bond
max_b_v(ii:785+1-ii) = params%maxi_bond

! Initializing b_v
DO ii=1,785
	IF (b_v(ii).GT.max_b_v(ii)) THEN
		b_v(ii) = max_b_v(ii)
	END IF
END DO

! Allocating vector of tensors
ALLOCATE(matrices%vts(784))

! Randomly initializing tensors
DO ii = 1, 784
	ALLOCATE(matrices%vts(ii)%elem( b_v(ii), 2, b_v(ii+1)) )
	! A tensors (using random values between 0 and 1)
	CALL RANDOM_NUMBER(matrices%vts(ii)%elem(:,:,:))
END DO

! Checking correct dimensions of bonds
CALL CHECK_BONDS(b_v, max_b_v)

! Applying LEFT_CANO subroutine
CALL LEFT_CANO(matrices%vts, b_v)

! Initializing cumulants
CALL INIT_CUMULANTS(dataset, matrices%vts, cumuls%vtm, b_v)


! ---------- GRADIENT DESCENT ----------
WRITE(6,*) 'GRADIENT_DESCENT'

! Saving parameters to file (conventionally doing it here)
OPEN(25, FILE="Saved/parameters.txt", ACTION="WRITE")
WRITE(25,*) params
CLOSE(25)

! Start interaging
DO ii=1, params%num_of_iter
	WRITE(*,*) "Step:", ii

	! Evaluating a right-to-left and a successive left-to-right GD
	CALL GRADIENT_DESCENT(b_v, cumuls%vtm, dataset, matrices%vts, params%learning_rate, &
			& params%cut_off, params%descent_steps, params%mini_bond, max_b_v)
	
	! Eventually make LR decrease during training to avoid
	!	jumping out from loss function minima
	params%learning_rate = params%learning_rate*0.9


	! Printing to terminal boundaries of tensors in a human readable way
	!	to see how the training procedure is going
	DO kk=1,28
		DO jj = 1,27
			WRITE(6,"(I3)", ADVANCE="No") b_v((kk-1)*28+jj)
		END DO
		WRITE(6,"(I3)", ADVANCE="Yes") b_v(kk*28)
	END DO
	
	! Saving bonds to file
	WRITE(num_file, "(I3.3)") ii
	OPEN(25, FILE="Saved/bonds_evol.txt", ACTION="WRITE", STATUS="Unknown", ACCESS="Append")
	DO kk=1,785
		WRITE(25,"(I3.3)", ADVANCE="No") b_v(kk)
		WRITE(25,"(A1)", ADVANCE="No") ","
	END DO
	WRITE(25,*) ""
	CLOSE(25)

	
	! SAVING PARAMETERS AT EVERY ITERATION
	WRITE(*,*) "Saving parameters and network..."

	! Saving bonds
	OPEN(25, FILE="Saved/bonds.txt", ACTION="WRITE")
	WRITE(25,*) b_v
	CLOSE(25)

	! Saving all tensors
	DO kk=1,784
		WRITE(num_file, "(I3.3)") kk
		OPEN(26, FILE="Saved/tensors/t"//num_file, ACTION="WRITE")
		WRITE(26,*) matrices%vts(kk)%elem
		CLOSE(26)
	END DO

	! Writing the number of epoch to file
	OPEN(26, FILE="Saved/Stop_epoch.txt", ACTION="WRITE")
	WRITE(26,"(I3.3)") ii
	CLOSE(26)
	
	WRITE(*,*) "Saved"
	
END DO

END PROGRAM ULTN_TRAIN