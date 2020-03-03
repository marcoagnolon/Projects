! PROGRAM TO TRAIN A Matrix Product States (MPS) OBJECT
!	ON A FRACTION OF MNIST DATASET (1000 SAMPLES) TO TRY 
!	TO LEARN PIXEL DISTRIBUTION FEATURES (UNSUPERVISED LEARNING TASK)

PROGRAM ULTN_GENERATE

! Importing useful modules
USE UTILITY

IMPLICIT NONE

! Defyning some variables
! dataset_load			: the set of images (70k * (28 * 28) array, with elements {0,1})
!							It is used to easily load the file, then train elements
!							are selected and used in "dataset" variable
! dataset				: the selected elements to train on
! dataset_mask_r1,2		: the indexes mask of random numbers
! dataset_mask_i1,2		: the indexes mask of integer numbers to select elements from dataset_load
!						N.B.: the masks are 2 to make sure that the second is different from the first
!								thus testing the MPS reconstruction power on new images
! b_v					: array representing the bonds values (dimension = 785, because it contains
!							the dimension of the left and right bonds of 784 middle tensors)
! max_b_v				: maximum bonds allowed (otherwise the tensors can not be reconstructed
!							after being merged)
! matrices				: list of tensors of the network
! img_reco				: variable used to store the generated/reconstructed images
!							then they are printed to files and the variable can be used again
! accuracy,acc			: variables to measure the accuracy
! ii, jj				: indexes to be used in loops
! chosen_seed			: array of number to use as seed to have reproducible executions
! num_file				: string useful to save several files

!-----DATASET VARIABLES-----
INTEGER(KIND=1), DIMENSION(70000,784) :: dataset_load
INTEGER(KIND=1), DIMENSION(2000,784) :: dataset
REAL*8, DIMENSION(1000) :: dataset_mask_r1
REAL*8, DIMENSION(2000) :: dataset_mask_r2
INTEGER, DIMENSION(1000) :: dataset_mask_i1
INTEGER, DIMENSION(2000) :: dataset_mask_i2
!-----MPS & GENERATE/RECONSTRUCTION VARIABLES-----
INTEGER, DIMENSION(785) :: b_v
TYPE(VEC_TENS) :: matrices
INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE :: img_reco
REAL*8 ::  accuracy, acc
!-----OTHER VARIABLES-----
INTEGER :: ii, jj
INTEGER, DIMENSION(33) :: chosen_seed
CHARACTER(LEN=3) :: num_file

! STARTING COMPUTATIONS

! Allocating vector of tensors
ALLOCATE(matrices%vts(784))

! Allocating img_reco
ALLOCATE(img_reco(784))
! Setting img_reco to 0
img_reco = 0

! Initializing seed using random but fixed values (priviously sampled)
chosen_seed = (/70290, 38879, 10013, 13932, 41407, 69718, 84713, 49520, 32977, 46194, 81813, 39405, &
& 85864, 19900, 73449, 56005, 19063, 74073, 25801, 63107, 71317, 22322, 90386, 40631, 99486, 65006, &
& 89833, 49354, 69803, 70049, 54845, 52086, 76196/)
CALL RANDOM_SEED(put=chosen_seed)

! Reading dataset from file
OPEN(33, FILE="../TRAIN_PROGRAM/data.bin",FORM="UNFORMATTED", STATUS="UNKNOWN", ACTION="READ", ACCESS='STREAM')
READ(33) dataset_load
WRITE(*,*) SHAPE(dataset)
CLOSE(33)

! Initializing 1000 random numbers in [0,1)
CALL RANDOM_NUMBER(dataset_mask_r1)
CALL RANDOM_NUMBER(dataset_mask_r2)

! Scaling: [0,1) → [0, 70k)
!	Truncating to integer
!	Adding 1 to all indexes
dataset_mask_i1 = 1+FLOOR(dataset_mask_r1*70000)
dataset_mask_i2 = 1+FLOOR(dataset_mask_r2*70000)

! Comparing masks to have all elements in the second mask different from the first mask
!	so that all the images to be reconstructed are different from those trained on
DO ii=1,SIZE(dataset_mask_i1)
	DO jj=1,SIZE(dataset_mask_i2)
		IF (dataset_mask_i1(ii).EQ.dataset_mask_i2(jj)) THEN
			! If there are match, changing the mask in this simple way
			dataset_mask_i2(jj) = dataset_mask_i2(jj) + 1
		END IF
	END DO
END DO

dataset = dataset_load(dataset_mask_i2,:)

! Verbose
WRITE(*,*) "Dataset is loaded."


! LOADING MPS OBJECT

! Loading dimensions
OPEN(25, FILE="ToLoad/bonds.txt", ACTION="READ")
READ(25,*) b_v
CLOSE(25)

! Allocating all tensors
DO ii = 1, 784
	ALLOCATE(matrices%vts(ii)%elem( b_v(ii), 2, b_v(ii+1)) )
END DO

! Initializing all tensors (reading them from file)
DO ii=1,784
	WRITE(num_file, "(I3.3)") ii
	OPEN(26, FILE="ToLoad/tensors/t"//num_file, ACTION="READ")
	READ(26,*) matrices%vts(ii)%elem
	CLOSE(26)
END DO

! Normalizing the last tensor to generate
matrices%vts(784)%elem(:,:,:) = matrices%vts(784)%elem(:,:,:)/SQRT(SUM(matrices%vts(784)%elem(:,:,:)**2))

! Verbose
WRITE(*,*) "Last tensor normalized:", SUM(matrices%vts(784)%elem(:,:,1)**2)

! Generating 10 images
DO jj=1,10

	! Generating image
	CALL GENERATE(b_v, matrices%vts, img_reco)

	! Converting index jj to string
	WRITE(num_file,"(I3.3)") jj

	! Saving img_reco to file
	OPEN(25, FILE="ToLoad/Img_gen_"//num_file//".txt", ACTION="WRITE")
	DO ii=1,783
		WRITE(25,"(I1.1)") img_reco(ii)
	END DO
	WRITE(25,"(I1.1)", ADVANCE="No") img_reco(784)
	CLOSE(25)

END DO


! Reconstructing 10 images
! Resetting accuracy
accuracy = 0.0

DO jj=1,2000

	! Verbose
	WRITE(*,*) "Img", jj

	! Reconstructing image
	CALL GENERATE_FROM_BOT(b_v, matrices%vts, dataset(jj,:), img_reco, acc)

	! Converting index jj to string
	!WRITE(num_file,"(I3.3)") jj

	! Saving img_reco to file
	!OPEN(25, FILE="ToLoad/Img_reco_"//num_file//".txt", ACTION="WRITE")
	!DO ii=1,783
	!	WRITE(25,"(I1.1)", ADVANCE="Yes") img_reco(ii)
	!END DO
	!WRITE(25,"(I1.1)", ADVANCE="No") img_reco(784)
	!CLOSE(25)
	
	! Updating accuracy
	accuracy = accuracy + acc

END DO

! Verbose
WRITE(*,*) 'Accuracy', accuracy/2000.0

! Writing accuracy to file
OPEN(25, FILE="ToLoad/Accuracy.txt", ACTION="WRITE")
WRITE(25,*) accuracy
CLOSE(25)

END PROGRAM ULTN_GENERATE