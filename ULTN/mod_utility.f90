! Defyning TYPES and FUNCTIONS to manipulate tensors
!	and to take out training, generating and reconstructing procedure
MODULE UTILITY
IMPLICIT NONE

	! TYPE to store computing parameters
	TYPE PARAMETERS
		! space_size			: numb. of input size (28x28)
		! descent_steps			: number of steps for gradient descent 
		! num_of_iter			: number of calls to gradient descent
		! init_bond_dim			: starting dimension of each bond
		! maxi_bond				: maximum bond dimension allowed
		! mini_bond				: minimum bond dimension allowed
		! cut_off				: threshold for SVD procedure
		! learning_rate			: gradient descent learning rate
		
		INTEGER :: space_size,  descent_steps, n_batch, num_of_iter, init_bond_dim, maxi_bond, mini_bond
		REAL*8 :: cut_off, learning_rate
		
	END TYPE PARAMETERS
	
	
	! TYPE for the ORDER 3 TENSOR OBJECT
	TYPE TENSOR3
		! elem				: middle tensor values
		REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: elem
	END TYPE TENSOR3

	
	! TYPE ARRAY OF TENSORS (trick to have RUN-TIME ALLOCATABLE tensors)
	TYPE VEC_TENS
		! vts				: vector of tensors
		TYPE(TENSOR3), DIMENSION(:), ALLOCATABLE :: vts
	END TYPE VEC_TENS


	! TYPE CUMULANTS
	!	Cumulants represent a trick to avoid repeated and identical tensor contractions
	TYPE CUMULANT_ELE
		! elem				: cumulant values
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: elem
	END TYPE CUMULANT_ELE
	
	
	! TYPE ARRAY OF CUMULANTS (trick to have RUN-TIME ALLOCATABLE objects)
	TYPE CUMULANTS
		! vtm				: vector of matrices
		TYPE(CUMULANT_ELE), DIMENSION(:), ALLOCATABLE :: vtm
	END TYPE CUMULANTS

! ------------------------------------------------------------------------------------------------------------

	! Defining interface to use MERGED_TENSOR function
	
	!INTERFACE MERGE_TENSOR
	!	MODULE PROCEDURE MERGED_TENSOR
	!END INTERFACE


	! DEFYNING FUNCTIONS
	CONTAINS

		! FUNCTION TO MERGE TWO TENSORS
		FUNCTION MERGED_TENSOR(t1, t2)
		!-----INPUT-----
			! t1, t2			: tensors to be merged
		!-----IN/OUT-----
			! ///
		!-----OUTPUT-----
			! MERGED_TENSOR		: tensor to be returned (the result of the merging)
		!-----VARIABLES-----
			! shape1,shape2		: vector of dimensions of t1%elem,t2%elem
			! ii,jj,kk,ll,mm	: indexes to be used in loops
			
			! ARGUMENTS
			TYPE(TENSOR3), INTENT(IN) :: t1, t2
			REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: MERGED_TENSOR
			! VARIABLES
			INTEGER, DIMENSION(3) :: shape1, shape2
			INTEGER :: ii,jj,kk,ll,mm
			
			! STARTING COMPUTATIONS
					
			shape1 = SHAPE(t1%elem)
			shape2 = SHAPE(t2%elem)
			
			ALLOCATE(MERGED_TENSOR(shape1(1), shape1(2), shape2(2),shape2(3) ))
			
			!IF (shape1(3) .NE. shape2(1))
			!	print*, 'dimensioality problem'
			!END IF

			! Initializing all elements to 0
			MERGED_TENSOR = 0
			
			DO ii = 1,shape1(1)
				DO jj = 1,shape1(2)
					DO ll = 1,shape2(2)
						DO mm = 1,shape2(3)
							DO kk = 1,shape1(3)
								MERGED_TENSOR(ii,jj,ll,mm)= MERGED_TENSOR(ii,jj,ll,mm) + t1%elem(ii,jj,kk)*t2%elem(kk,ll,mm)
							END DO
						END DO
					END DO
				END DO
			END DO
			
			RETURN
			
		END FUNCTION MERGED_TENSOR



		! SUBROUTINE TO SPLIT A MERGED TENSOR INTO 2 TENSORS (DEPENDING ON THE GIVEN THRESHOLD or ON keep_dim ARGUMENT)
		SUBROUTINE SPLIT_TENSOR(tensor_merged, U_out, VT_out, THRESH, going_right, mini_dim, maxi_dim, previous_bond, keep_dim)
		!-----INPUT-----
			! tensor_merged			: tensor to be split
			! THRESH				: fraction of the largest singolar value (S[0]) used to decide which other singolar values to keep:
			!							s(ii) is kept if s(ii) > s(0)*THRESH
			! going_right			: boolean variable that states if the process (gradient descent in our case) is
			!							going right in the array of tensors
			!							- When going right, normalized eigenvectors are on the left
			!								(eigenvalues are multipied on right)
			!							- when going left, normalized eigenvectors is on the right
			! mini_dim				: value representing the minimum dimension (used to start the cut procedure)
			! maxi_dim				: value representing the maximum dimension of the bonds
			! keep_dim				: boolean that states if act the cut using the threshold or not
			!							(not acting the cut is useful in the generative process)
		!-----IN/OUT-----
			! U_out					: it is the left tensor; it is modified on output
			! V_out					: it is the right tensor; it is modified on output
			! previous_bond			: it is the dimension (corresponding to "k" index in the notation A_ijk · B_klm)
			!							of the tensors which had been contracted
		!-----OUTPUT-----
			! ///
		!-----VARIABLES-----
			! tensor_merged_W		: tensor to be used to make the evaluations (it has different shape w.r.t. INPUT tensor_merged)
			! ii					: index to be used in loops
			! cont_val				: number of singular values to keep (depending on evaluations based on threshold)
			! size_S				: size of the array of the singular values resulting from SVD procedure
			! M						: vector of dimensions of tensor_merged
			! U_trunc,VT_trunc,U,VT	: variables useful during progressive evaluations steps
			
			! OTHERS				: variables useful to make llapack subroutine act properly
			
			! ARGUMENTS
			REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(IN) :: tensor_merged
			REAL*8, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: U_out
			REAL*8, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: VT_out
			REAL*8, INTENT(IN) :: THRESH
			LOGICAL, INTENT(IN) :: going_right, keep_dim
			INTEGER, INTENT(IN) :: mini_dim, maxi_dim
			INTEGER, INTENT(INOUT) :: previous_bond
			
			! VARIABLES
			REAL*8, DIMENSION(:,:), ALLOCATABLE :: tensor_merged_W
			INTEGER :: cont_val, size_S, ii
			INTEGER, DIMENSION(SIZE(SHAPE(tensor_merged))) :: M
			REAL*8, DIMENSION(:,:), ALLOCATABLE :: U_trunc, VT_trunc
			REAL*8, DIMENSION(:,:), ALLOCATABLE :: U, VT
			
			! VARIABLES USEFUL FOR LLAPACK SUBROUTINE
			REAL*8, DIMENSION(:), ALLOCATABLE:: WORK, S
			INTEGER :: LDA, LDU, LDVT, LWORK, INFO
			INTEGER, DIMENSION(:), ALLOCATABLE:: IWORK
			
			! STARTING COMPUTATIONS
			
			! N.B.: for our purposes, it should be M(2) = M(3) = 2 (because they are input dimensions)
			M=SHAPE(tensor_merged)
			
			ALLOCATE(tensor_merged_W(M(1)*M(2),M(3)*M(4)))
			tensor_merged_W = RESHAPE(tensor_merged, (/M(1)*M(2),M(3)*M(4)/))

					
			!-----COMPUTING EIGENVALUES AND EIGENVECTORS-----
			
			ALLOCATE(U(M(1)*M(2),M(1)*M(2)))
			ALLOCATE(VT(M(3)*M(4),M(3)*M(4)))
			
			! Link to DGESDD docum:
			!	http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_gad8e0f1c83a78d3d4858eaaa88a1c5ab1.html#gad8e0f1c83a78d3d4858eaaa88a1c5ab1
			! Initializing variables following the above documentation
			LDA = M(1)*M(2)
			LDU = M(1)*M(2)
			LDVT = M(3)*M(4)
			ALLOCATE(S(MIN(M(1)*M(2), M(3)*M(4))))
			ALLOCATE(WORK(1))
			ALLOCATE(IWORK( MAX(1,8*MIN(M(1)*M(2),M(3)*M(4))) ) )

			! Calling subroutines to initialize correct dimensions
			LWORK=-1
			CALL DGESDD('A', M(1)*M(2), M(3)*M(4), tensor_merged_W, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)

			! Initilialize correct dimensions
			LWORK = INT(WORK(1))
			DEALLOCATE(WORK)
			ALLOCATE(WORK(LWORK))
			! Recalling subroutine with correct dimensions
			CALL DGESDD('A', M(1)*M(2), M(3)*M(4), tensor_merged_W, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO)
			
			! Retrieving S size
			size_S = SIZE(S)
			
			! Evaluating cont_val
			IF (keep_dim) THEN
				! Keeping the maximum number of indexes allowed by the system
				cont_val = MIN(previous_bond,size_S)
				! It should be
				!		cont_val = previous_bond
				! but sometimes it can happen that the nearer bonds do not allow to retrieve previous_bond
			ELSE
				! Counting how many singolar values to keep using THRESH
				cont_val = 1
				DO WHILE (cont_val.LE.maxi_dim .AND. cont_val.LE.size_S)
					IF (S(cont_val).GE.S(1)*THRESH) THEN
						cont_val = cont_val + 1
					ELSE
						! If there are no valid singular values to keep, then the loop can exit
						EXIT
					END IF
				END  DO
				! Fixing number of singolar values (mismatch between index and counter
				!	due to the fact that cont_val is incremented even after the last test)
				cont_val = MAX(mini_dim, cont_val-1)
				! Updating the bond value
				previous_bond = cont_val
			END IF
			
			ALLOCATE(U_trunc(M(1)*M(2), cont_val))
			ALLOCATE(VT_trunc(cont_val,M(3)*M(4)))

			! Evaluating output variables (depending on going_right variable)
			IF (going_right) THEN
				! Truncating U
				U_trunc(:,:) = U(:, :cont_val)
				! Truncating and multiplying rows of VT
				DO ii=1, cont_val
					VT_trunc(ii,:) = S(ii) * VT(ii, :)
				END DO
				! "Normalizing" VT_trunc
				VT_trunc = VT_trunc / DSQRT(SUM(VT_trunc**2))
			ELSE
				! Truncating and multiplying colums of U
				DO ii=1, cont_val
					U_trunc(:,ii) = U(:, ii) * S(ii)
				END DO
				! "Normalizing" U_trunc
				U_trunc = U_trunc / DSQRT(SUM(U_trunc**2))
				! Truncating VT
				VT_trunc(:,:) = VT(:cont_val, :)
			END IF
					
			! ALLOCATING OUTPUT TENSORS (first deleting them)
			DEALLOCATE(U_out)
			DEALLOCATE(VT_out)
			ALLOCATE(U_out(M(1), M(2), cont_val))
			ALLOCATE(VT_out(cont_val, M(3), M(4)))
			
			! Reshaping output tensors
			U_out = RESHAPE(U_trunc, (/M(1), M(2), cont_val/))
			VT_out = RESHAPE(VT_trunc, (/cont_val, M(3), M(4)/))

		END SUBROUTINE SPLIT_TENSOR



		! SUBROUTINE TO MAKE THE TENSORS ARRAY LEFT CANONICALIZED (rightmost tensor has not
		!		normalized eigenvectors, but it is normalized w.r.t. its norm)
		SUBROUTINE LEFT_CANO(tensors_list, bonds_list)
		!-----INPUT-----
			! ///
		!-----IN/OUT-----
			! tensors_list		: the array of MPS tensors to be left-canonicalized
			! bonds_list		: the array of MPS tensors boundaries to be used for tensors reconstructions
		!-----OUTPUT-----
			! ///
		!-----VARIABLES-----
			! shape1, shape2	: the shapes of the tensors to be split
		
			! ARGUMENTS
			TYPE(TENSOR3), DIMENSION(:), INTENT(INOUT) :: tensors_list
			INTEGER, DIMENSION(:), INTENT(INOUT) :: bonds_list
			
			! VARIABLES
			INTEGER, DIMENSION(3) :: shape1, shape2

			! STARTING COMPUTATIONS

			INTEGER :: ii
			
			DO ii=1,SIZE(bonds_list)-2		! bonds_list has one more bond (that at the end)
				
				! Retrieving shapes
				shape1 = SHAPE(tensors_list(ii)%elem)
				shape2 = SHAPE(tensors_list(ii+1)%elem)
				
				! Conventionally setting THRESH=1, mini_bond=1, maxi_bond=1 because when "keep_dim"=.TRUE. they are ignored
				CALL SPLIT_TENSOR( MERGED_TENSOR(tensors_list(ii), tensors_list(ii+1)), &
			&		tensors_list(ii)%elem, tensors_list(ii+1)%elem, 0d0, .TRUE., 1, 1, bonds_list(ii+1), .TRUE.)
				
			END DO
			
		END SUBROUTINE LEFT_CANO
	


		! SUBROUTINE TO CHECK THAT THE CURRENT BONDS ARE ALL ≤ max_bonds
		SUBROUTINE CHECK_BONDS(curr_bonds, maxi_bonds)
		!-----INPUT-----
			! curr_bonds	: current bonds (these values are modified during training)
			! maxi_bonds	: these are the maximum values allowed and choosen by the user
		!-----IN/OUT-----
			! ///
		!-----OUTPUT-----
			! ///
		!-----VARIABLES-----
			! ii			: index to be used in loops
			
			! ARGUMENTS
			INTEGER, DIMENSION(:), INTENT(IN) :: curr_bonds, maxi_bonds
			! VARIABLES
			INTEGER :: ii
			
			! STARTING COMPUTATIONS
			
			IF (SIZE(curr_bonds).NE.SIZE(maxi_bonds)) THEN
				WRITE(6,*) "ERROR: the dimensions of the CURRENT BONDS and MAXIMUM BONDS vectors do not match..."
				WRITE(6,*) "Ending the process..."
				STOP
			END IF
			
			DO ii=1,SIZE(curr_bonds)
				IF (curr_bonds(ii).GT.maxi_bonds(ii)) THEN
					WRITE(6,*) "ERROR: some value of the CURRENT BONDS are LARGER than the corresponding MAXIMUM BONDS."
					WRITE(6,*) "Ending the process..."
					STOP
				END IF
			END DO
		
		! IF THIS SUBROUTINE RUNS WITH NO ERRORS, THEN THE CONDITIONS ARE RESPECTED
		
		END SUBROUTINE CHECK_BONDS


		! SUBROUTINE TO INITIALIZE THE CUMULANTS
		!	N.B.: SUPPOSING THAT THE TRAIN STARTS FROM THE RIGHTMOST TENSOR
		SUBROUTINE INIT_CUMULANTS(dataset, tensors, cumuls, bonds_list)
		!-----INPUT-----
			! dataset			: the array of input images; they are used to evaluate the cumulants
			! tensors			: the array of MPS tensors; they are used to evaluate the cumulants
			! bonds_list		: the array of MPS tensors boundaries
		!-----IN/OUT-----
			! ///
		!-----OUTPUT-----
			! cumuls			: the initialized cumulants
		!-----VARIABLES-----
			! jj, kk, nn		: indexes to be used in loops
			! Nx				: the number of input images (i.e.: the size of the dataset)
			
			! ARGUMENTS
			INTEGER(KIND=1), DIMENSION(:,:), INTENT(IN) :: dataset
			TYPE(TENSOR3), DIMENSION(:), INTENT(IN) :: tensors
			TYPE(CUMULANT_ELE), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: cumuls
			INTEGER, DIMENSION(:), INTENT(IN) :: bonds_list
			
			! VARIABLES
			INTEGER :: jj,kk,nn
			INTEGER, DIMENSION(2) :: Nx

			! STARTING COMPUTATIONS

			Nx=SHAPE(dataset)

			! Allocating 784 cumulants
			ALLOCATE(cumuls(784))

			! Allocating each cumulant with proper dimensions
			!	N.B.: the process is supposed to start from right
			!			so the rightmost cumulant will be transposed
			DO nn=1,783
				ALLOCATE(cumuls(nn)%elem(Nx(1), bonds_list(nn)) )
				cumuls(nn)%elem = 0.0
			END DO
			
			! Allocating the rightmost cumulant (with transposed dimensions
			ALLOCATE(cumuls(784)%elem(bonds_list(785), Nx(1)) )
			
			! The first and the last cumulants are "null elements"
			!	thus they are arrays of 1s.
			cumuls(1)%elem(:,:) = 1.0
			cumuls(784)%elem(:,:) = 1.0
			
			! Evaluating the cumulants starting from left
			DO nn = 1, 782
				DO kk = 1, bonds_list(nn+1)
					DO jj = 1, bonds_list(nn)
						cumuls(nn+1)%elem(:,kk) = cumuls(nn+1)%elem(:,kk) + cumuls(nn)%elem(:,jj) &
					&		 * tensors(nn)%elem(jj,dataset(:,nn)+1,kk)
								!								^^ Here the "+1" is to have exploit the fact that
								!						the dataset contains values {0,1}, so to use them as indexes
					END DO
				END DO
			END DO
			
		END SUBROUTINE INIT_CUMULANTS



		! SUBROUTINE TO UPDATE CUMULANTS
		SUBROUTINE UPDATE_CUMULANTS(up_cumuls, up_tensors, bonds_list_up, dataset, c_idx, going_right_upd)
		!-----INPUT-----
			! up_tensors			: the array of MPS tensors (they are not be modified by this subroutine)
			! dataset				: the array of input images
		!-----IN/OUT-----
			! bonds_list_up			: the array of MPS tensors boundaries
			! up_cumuls				: cumulants (in this subroutine they have a different name)
			! c_idx					: the current index corresponding to the merged bond
			! going_right_upd		: boolean that states if the training procedure is going right,
			!							the update happens in different ways basing on this variable
		!-----OUTPUT-----
			! ///
		!-----VARIABLES-----
			! Nx					: the number of input images (i.e.: the size of the dataset)
			! jj, kk				: indexes to be used in loops
		
			TYPE(CUMULANT_ELE), DIMENSION(:), INTENT(INOUT) :: up_cumuls
			TYPE(TENSOR3), DIMENSION(:), INTENT(IN) :: up_tensors
			INTEGER, DIMENSION(:), INTENT(INOUT) :: bonds_list_up
			INTEGER(KIND=1),  DIMENSION(:,:), INTENT(IN) :: dataset
			INTEGER, INTENT(IN) :: c_idx
			LOGICAL, INTENT(IN) :: going_right_upd
			
			INTEGER, DIMENSION(2) :: Nx
			INTEGER :: jj, kk
			
			! STARTING COMPUTATIONS
			
			Nx=SHAPE(dataset)

			! GOING RIGHT
			IF (going_right_upd) THEN
				! Deleting cumulant on the left (index = c_idx)
				DEALLOCATE(up_cumuls(c_idx)%elem)
				! Allocating a cumulant with "left shape": [ N_samples, bond(current indx) ]
				ALLOCATE(up_cumuls(c_idx)%elem(Nx(1), bonds_list_up(c_idx)))				
				up_cumuls(c_idx)%elem = 0
				
				! Computing new cumulant using the left ones
				DO kk = 1, bonds_list_up(c_idx)
					DO jj = 1, bonds_list_up(c_idx-1)
						! cumul(curr_idx) = cumul(curr_idx-1) X A(curr_idx-1)
						up_cumuls(c_idx)%elem(:,kk) = up_cumuls(c_idx)%elem(:,kk) + up_cumuls(c_idx-1)%elem(:,jj) &
					&		 * up_tensors(c_idx-1)%elem(jj,dataset(:,c_idx-1)+1,kk)
								!											 ^^ Here the "+1" is to have exploit the fact that
								!								the dataset contains values {0,1}, so to use them as indexes
					END DO
				END DO

			! GOING LEFT
			ELSE
				! Deleting cumulant on the left (index = c_idx-1)
				DEALLOCATE(up_cumuls(c_idx-1)%elem)
				! Allocating a cumulant with "right shape": [ bond(current indx), N_samples ]
				ALLOCATE(up_cumuls(c_idx-1)%elem(bonds_list_up(c_idx), Nx(1)))
				up_cumuls(c_idx-1)%elem = 0

				! Computing new cumulant using the right ones
				DO jj = 1, bonds_list_up(c_idx)
					DO kk = 1, bonds_list_up(c_idx+1)
						! cumul^T(curr_idx-1) = cumul(curr_idx) X A(curr_idx)
						up_cumuls(c_idx-1)%elem(jj,:) = up_cumuls(c_idx-1)%elem(jj,:) +&
					&		up_tensors(c_idx)%elem(jj,dataset(:,c_idx)+1,kk) * up_cumuls(c_idx)%elem(kk,:)
					END DO
				END DO
			END IF
			
		END SUBROUTINE UPDATE_CUMULANTS


! -----------------------------------------------------------------------------------------------------------------------------------------
		! SUBROUTINE THAT UPDATES THE MPS TENSORS USING GRADIENT DESCENT
		! 	Starting from the rightmost pair of MPS tensors, it MERGES,UPDATES,SPLITS
		!		them iteratively, going left and then returning right
		SUBROUTINE GRADIENT_DESCENT(bonds_list, cumuls, dataset, tensors, learning_rate, split_tresh, &
		& steps_num, mini_dim_grad, maxi_dim_grad)
		!-----INPUT-----
			! dataset				: the array of input images
			! learning_rate			: the step lenght to be done at every GD iteration
			! split_tresh			: the threshold passed when calling SPLIT_TENSOR subroutine
			! steps_num				: the number of iteration used for the update before split
			! mini_dim_grad			: the minimum dimension to pass to SPLIT_TENSOR subroutine
			! maxi_dim_grad			: the maximum dimension to pass to SPLIT_TENSOR subroutine
		!-----IN/OUT-----
			! bonds_list			: the list of the MPS tensors boundaries
			! cumuls				: the array of cumulants
			! tensors				: the array of MPS tensors
		!-----OUTPUT-----
			! ///
		!-----VARIABLES-----
			! ii, jj, kk			: indexes to be used in loops
			! pixel					: index representing the left tensors in the pairs of
			!							tensors when merging them to apply SGD
			! bond_idx (=pixel+1)	: index representing the bond around which the merging
			!							and splitting are applied
			! gd					: counter for GD loops
			! merged				: tensor resulting from the MERGED_TENSOR function applied
			!							to a pair of MPS tensor
			! gradient				: the tensor used to update the merged tensor
			! Nx					: the number of input images (i.e.: the size of the dataset)
			! dev_psi				: tensor representing the derivative of psi w.r.t. the merged tensor
			! psi,psi_inv			: arrays representing the wave function of the system and its inverse
			! left_vecs, right_vecs	: tensors representing the nearest cumulants (on left and on right)
			! loss					: the log-likelihood for the current psi  (evaluated at the end of the GD)
			! app
			
			! ARGUMENTS
			INTEGER, DIMENSION(:), INTENT(INOUT) :: bonds_list
			TYPE(CUMULANT_ELE), DIMENSION(:), INTENT(INOUT) :: cumuls
			INTEGER(KIND=1), DIMENSION(:,:), INTENT(IN) :: dataset
			TYPE(TENSOR3), DIMENSION(:), INTENT(INOUT) :: tensors
			REAL*8, INTENT(IN) :: learning_rate, split_tresh
			INTEGER, INTENT(IN) :: steps_num, mini_dim_grad
			INTEGER, DIMENSION(:), INTENT(IN) :: maxi_dim_grad
			
			! VARIABLES
			INTEGER :: ii, jj, kk, pixel, bond_idx, gd
			REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: merged, gradient
			INTEGER, DIMENSION(2) :: Nx
			REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: dev_psi
			REAL*8, DIMENSION(:), ALLOCATABLE :: psi, psi_inv
			REAL*8, DIMENSION(:,:), ALLOCATABLE :: left_vecs, right_vecs
			REAL*8 :: loss, app
			app = 0d0
			
			! STARTING COMPUTATIONS
			
			Nx=SHAPE(dataset)
			

			ALLOCATE(psi(Nx(1)))
			ALLOCATE(psi_inv(Nx(1)))

			! Resetting loss
			loss = 0d0
			
			! GD GOING LEFT
			DO pixel = 783, 1, -1
				
				! retrieving the bond index
				bond_idx=pixel+1

				! evaluating the merged tensor
				merged=MERGED_TENSOR(tensors(pixel), tensors(pixel+1))

				! Repeating GD
				DO gd=1, steps_num
					ALLOCATE(left_vecs(Nx(1),bonds_list(bond_idx-1)))
					ALLOCATE(right_vecs(bonds_list(bond_idx+1),Nx(1)))			
					ALLOCATE(dev_psi(Nx(1),bonds_list(bond_idx-1),bonds_list(bond_idx+1)))
					
					! (Re)Setting all to 0s
					dev_psi = 0d0
					left_vecs = 0d0
					right_vecs = 0d0
					psi = 0d0
					psi_inv = 0d0

					! Retrieving cumulants
					left_vecs = cumuls(bond_idx-1)%elem
					right_vecs = cumuls(bond_idx)%elem

					! Evaluating psi and its derivative w.r.t. the merged tensor
					DO ii=1, Nx(1)
						DO jj=1, bonds_list(bond_idx-1)
							DO kk=1, bonds_list(bond_idx+1)
								dev_psi(ii,jj,kk)=dev_psi(ii,jj,kk) + left_vecs(ii,jj)*right_vecs(kk,ii)
								psi(ii) =  psi(ii) + left_vecs(ii,jj) * merged(jj, dataset(ii,pixel)+1, &
							&				dataset(ii, pixel+1)+1, kk) * right_vecs(kk,ii)
							END DO
						END DO
					END DO

					! Evaluating the inverse of psi
					psi_inv = 1d0/psi

					! (Re)setting gradient
					ALLOCATE(gradient(bonds_list(bond_idx-1),2,2,bonds_list(bond_idx+1)))
					gradient = 0d0
					
					! Evaluating gradient
					DO ii=1, 2
						DO jj=1, 2
							DO kk=1,Nx(1)
								IF ((dataset(kk,pixel).EQ.(ii-1)) .AND. (dataset(kk,pixel+1).EQ.(jj-1))) THEN
									gradient(:,ii,jj,:) = gradient(:,ii,jj,:) + dev_psi(kk,:,:) * psi_inv(kk)
								END IF
							END DO
						END DO
					END DO

					! Completing gradient evaluation
					! N.B.: thanks to the left_cano subroutine, here there is no need
					!		to divide by Z (the partition function)
					gradient = 2d0*gradient / Nx(1) - 2d0 * merged
					
					! Updating "merged" tensor
					merged = merged + gradient*learning_rate
										
					! Re-normalizing the merged tensor
					merged = merged / DSQRT(SUM(merged**2))

					! Deallocating
					DEALLOCATE(left_vecs)
					DEALLOCATE(right_vecs)
					DEALLOCATE(dev_psi)
					
					DEALLOCATE(gradient)
					
				END DO

				! Splitting the updated merged tensor
				CALL SPLIT_TENSOR(merged, tensors(pixel)%elem, tensors(pixel+1)%elem, split_tresh, .FALSE.,&
			&							mini_dim_grad, maxi_dim_grad(pixel+1), bonds_list(pixel+1), .FALSE.)
				DEALLOCATE(merged)
				
				! Updating cumulants
				IF (bond_idx.GT.2) THEN		! Can not update the first cumulant (it has to remain array of 1s)
					CALL UPDATE_CUMULANTS(cumuls, tensors, bonds_list, dataset, bond_idx, .FALSE.)
				END IF
				
			END DO

			! Printing boundaries of tensors to see how the training procedure is going
			!	(also formatting them in a 28*28 fancy way to print)
			DO kk=1,28
				DO jj = 1,27
					WRITE(6,"(I3)", ADVANCE="No") bonds_list((kk-1)*28+jj)
				END DO
				WRITE(6,"(I3)", ADVANCE="Yes") bonds_list(kk*28)
			END DO
			
			! Verbose
			WRITE(6,*) "Half"

			! Printing current likelihood
			WRITE(6,*) "Loss", 2*SUM(-DLOG(ABS(psi)))/Nx(1)


			! GD GOING RIGHT
			DO pixel = 1, 783
				
				bond_idx=pixel+1
				
				merged=MERGED_TENSOR(tensors(pixel), tensors(pixel+1))

				! Repeating GD
				DO gd=1, steps_num
				
					ALLOCATE(left_vecs(Nx(1),bonds_list(bond_idx-1)))
					ALLOCATE(right_vecs(bonds_list(bond_idx+1),Nx(1)))			
					ALLOCATE(dev_psi(Nx(1),bonds_list(bond_idx-1),bonds_list(bond_idx+1)))
					
					! (Re)Setting all to 0s
					dev_psi = 0d0
					left_vecs = 0d0
					right_vecs = 0d0
					psi = 0d0
					psi_inv = 0d0
					
					! Retrieving cumulants
					left_vecs = cumuls(pixel)%elem
					right_vecs = cumuls(pixel+1)%elem   

					! Evaluating psi and its derivative w.r.t. the merged tensor
					DO ii=1, Nx(1)
						DO jj=1, bonds_list(bond_idx-1)
							DO kk=1, bonds_list(bond_idx+1)
								dev_psi(ii,jj,kk)=dev_psi(ii,jj,kk) + left_vecs(ii,jj)*right_vecs(kk,ii)
								psi(ii) =  psi(ii) + left_vecs(ii,jj) * merged(jj, dataset(ii,pixel)+1, dataset(ii, pixel+1)+1, kk) * right_vecs(kk,ii)
							END DO
						END DO
					END DO
					
					! Evaluating the inverse of psi
					psi_inv = 1d0/psi

					! (Re)setting gradient
					ALLOCATE(gradient(bonds_list(bond_idx-1),2,2,bonds_list(bond_idx+1)))
					gradient = 0d0
					
					! Evaluating gradient
					DO ii=1, 2
						DO jj=1, 2
							DO kk=1,Nx(1)
								IF ((dataset(kk,pixel).EQ.(ii-1)) .AND. (dataset(kk,pixel+1).EQ.(jj-1))) THEN
									gradient(:,ii,jj,:) = gradient(:,ii,jj,:) + dev_psi(kk,:,:) * psi_inv(kk)
								END IF
							END DO
						END DO
					END DO

					! Completing gradient evaluation
					! N.B.: thanks to the left_cano subroutine, here there is no need
					!			to divide by Z (the partition function)
					gradient = 2d0*gradient / Nx(1) - 2d0 * merged

					! Updating "merged" tensor
					merged = merged + gradient*learning_rate
										
					! Re-normalizing the merged tensor
					merged = merged / DSQRT(SUM(merged**2))
					
					! Deallocating
					DEALLOCATE(left_vecs)
					DEALLOCATE(right_vecs)
					DEALLOCATE(dev_psi)
					
					DEALLOCATE(gradient)
					
				END DO
				
				! Splitting the updated merged tensor
				CALL SPLIT_TENSOR(merged, tensors(pixel)%elem, tensors(pixel+1)%elem, split_tresh, .TRUE.,&
			&							mini_dim_grad, maxi_dim_grad(pixel+1), bonds_list(pixel+1), .FALSE.)
				DEALLOCATE(merged)
				
				! Updating cumulants
				IF (bond_idx.LT.784) THEN
					CALL UPDATE_CUMULANTS(cumuls, tensors, bonds_list, dataset, bond_idx, .TRUE.)
				END IF

			END DO

			! Evaluating loss (on the last psi)
			loss = 2*SUM(-DLOG(ABS(psi)))/Nx(1)

			! Printing current likelihood
			WRITE(*,*) "Current loss", loss
			
			! Saving loss into a file
			OPEN(33, FILE="Saved/loss.txt", ACTION='WRITE', STATUS='UNKNOWN', ACCESS = 'append')
			WRITE(33,*) loss
			CLOSE(33)
			
			! Resetting loss for next iteration
			loss=0d0

			! Deallocating
			DEALLOCATE(psi)
			DEALLOCATE(psi_inv)
		
		END SUBROUTINE GRADIENT_DESCENT



		! SUBROUTINE TO GENERATE FROM SCRATCH A CONFIGURATION
		SUBROUTINE GENERATE(bonds_list, tensors, img)
		!-----INPUT-----
			! bonds_list		: the array of MPS tensors boundaries
			! tensors			: the array of MPS tensors
		!-----IN/OUT-----
			! ///
		!-----OUTPUT-----
			! img				: the generated image
		!-----VARIABLES-----
			! ii, jj, kk		: indexes to be used in loops
			! prob				: variable used to sample numbers in [0,1)
			! vec				: vector representing the state of previous pixel
			! vec_act			: vector representing the state of current pixel
			!						N.B.: the norm of vec is restored to 1 in every iteration.
			!							The ratio between the norm of the 2 vectors
			!							represents the probability for a given pixel
			!							to be 0 or 1 (the standard test is made for 0,
			!							eventually calcoulos are repeated for the other case)
			
			! ARGUMENTS
			INTEGER, DIMENSION(:), INTENT(IN) :: bonds_list
			TYPE(TENSOR3), DIMENSION(:), INTENT(IN) :: tensors
			INTEGER(KIND=1),  DIMENSION(784), INTENT(OUT) :: img
			
			! VARIABLES
			INTEGER :: ii, jj, kk
			REAL*8 :: prob
			REAL*8, DIMENSION(:), ALLOCATABLE :: vec, vec_act
			
			! STARTING COMPUTATIONS
			
			! Sampling a value in [0,1)
			CALL RANDOM_NUMBER(prob)
			
			! vec should start as a vector of a single "1",
			!	thus multiplying by it the first time gives no effect
			
			! Chosing and fixing the first bit
			IF (prob.LE.DSQRT(SUM(tensors(784)%elem(:,1,:)**2)/SUM(tensors(784)%elem(:,:,:)**2))) THEN
				img(784) = 0
			ELSE
				img(784) = 1
			END IF
			
			! ALLOCATE vec
			ALLOCATE(vec(bonds_list(784)))
			
			! Sampling all bit one by one in successive steps
			DO ii=783,1,-1
			
				! Allocating vec_act
				ALLOCATE(vec_act(bonds_list(ii)))
				
				! Sampling a value in [0,1)
				CALL RANDOM_NUMBER(prob)
				
				! Resetting new vector
				vec_act = 0
				
				! Creating new vector
				DO jj=1,bonds_list(ii)
					DO kk=1,bonds_list(ii+1)
						vec_act(jj) = vec_act(jj) + tensors(ii)%elem(jj,1,kk) * vec(kk)
					END DO
				END DO
				
				! Verbose
				WRITE(*,*) "Vec_act norm (=Prob[pixel=0]): ", SUM(vec_act**2)
				
				! Chosing and fixing the current bit
				IF (prob.LE.SUM(vec_act**2)/SUM(vec**2)) THEN
					img(ii) = 0
					! Saving new vector
					DEALLOCATE(vec)
					ALLOCATE(vec(SIZE(vec_act)))
					vec = vec_act
				ELSE
					img(ii) = 1
					! Re-creating new vector
					vec_act = 0
					DO jj=1,bonds_list(ii)
						DO kk=1,bonds_list(ii+1)
							vec_act(jj) = vec_act(jj) + tensors(ii)%elem(jj,2,kk) * vec(kk)
						END DO
					END DO
					! Saving new vector
					DEALLOCATE(vec)
					ALLOCATE(vec(SIZE(vec_act)))
					vec = vec_act
				END IF

				! Re-normalizing VEC so that the numbers do not become too small
				vec = vec / DSQRT(SUM(vec**2))

				! Dellocating vec_act
				DEALLOCATE(vec_act)

			END DO
		
		END SUBROUTINE GENERATE


		! SUBROUTINE TO GENERATE A CONFIGURATION FROM THE BOTTOM HALF OF THE IMAGE
		SUBROUTINE GENERATE_FROM_BOT(bonds_list, tensors, bot_half, img, accuracy)
		!-----INPUT-----
			! bonds_list		: the array of MPS tensors boundaries
			! tensors			: the array of MPS tensors
			! bot_half			: the input image to be reconstructed (only its
			!						bottom half will be used in this subroutine)
		!-----IN/OUT-----
			! ///
		!-----OUTPUT-----
			! img				: the reconstructed image
			! accuracy			: the pixel averaged correct guess (/392 guessed)
		!-----VARIABLES-----
			! ii, jj, kk		: indexes to be used in loops
			! choice			: flag that represents which state to select (0 or 1)
			!						during the image copy
			! prob				: variable used to sample numbers in [0,1)
			! vec, vec_act		: see SUBROUTINE GENERATE
			
			! ARGUMENTS
			INTEGER, DIMENSION(:), INTENT(IN) :: bonds_list
			TYPE(TENSOR3), DIMENSION(:), INTENT(IN) :: tensors
			INTEGER(KIND=1),  DIMENSION(784), INTENT(IN) :: bot_half
			INTEGER(KIND=1),  DIMENSION(784), INTENT(OUT) :: img
			REAL*8 , INTENT(OUT) :: accuracy
			
			! VARIABLES			
			INTEGER :: ii, jj, kk
			INTEGER :: choice
			REAL*8 :: prob
			REAL*8, DIMENSION(:), ALLOCATABLE :: vec, vec_act
			
			! STARTING COMPUTATIONS
			
			! Sampling a value in [0,1)
			CALL RANDOM_NUMBER(prob)
			
			! ALLOCATE vec
			ALLOCATE(vec(bonds_list(784)))
			
			! FILLING (COPYING) THE BOTTOM HALF OF THE IMAGE
			DO ii=784,393,-1
				img(ii) = bot_half(ii)

				! Evaluating progressive probabilities
				vec = 1
				! Allocating vec_act
				ALLOCATE(vec_act(bonds_list(ii)))
				vec_act = 0
				! Building next vector (no need to extract probabilities,
				! 		the values are clamped)
				choice = bot_half(ii)+1

				! Creating new vector
				DO jj=1,bonds_list(ii)
					DO kk=1,bonds_list(ii+1)
						vec_act(jj) = vec_act(jj) + tensors(ii)%elem(jj,choice,kk) * vec(kk)
					END DO
				END DO
				! Saving new vector
				DEALLOCATE(vec)
				ALLOCATE(vec(SIZE(vec_act)))
				vec = vec_act
				! Re-normalizing VEC so that the numbers do not become too small
				vec = vec / DSQRT(SUM(vec**2))
				DEALLOCATE(vec_act)				
			END DO
			
			! Here img is half initialized, vec is the vec at pixel numb. 393
			
			! Repeating same procedure as in SUBROUTINE GENERATE
			DO ii=392,1,-1
			
				! Allocating vec_act
				ALLOCATE(vec_act(bonds_list(ii)))

				! Sampling a value in [0,1)
				CALL RANDOM_NUMBER(prob)

				! Resetting "new" vector
				vec_act = 0

				! Creating new vector
				DO jj=1,bonds_list(ii)
					DO kk=1,bonds_list(ii+1)
						vec_act(jj) = vec_act(jj) + tensors(ii)%elem(jj,1,kk) * vec(kk)
					END DO
				END DO

				! Verbose
				!WRITE(*,*) "Vec_act norm (=Prob[pixel=0]): ", SUM(vec_act**2)
				
				! Chosing and fixing the current bit
				IF (prob.LE.SUM(vec_act**2)/SUM(vec**2)) THEN
					img(ii) = 0
					! Saving new vector
					DEALLOCATE(vec)
					ALLOCATE(vec(SIZE(vec_act)))
					vec = vec_act
				ELSE
					img(ii) = 1
					! Re-creating new vector
					vec_act = 0
					DO jj=1,bonds_list(ii)
						DO kk=1,bonds_list(ii+1)
							vec_act(jj) = vec_act(jj) + tensors(ii)%elem(jj,2,kk) * vec(kk)
						END DO
					END DO
					! Saving new vector
					DEALLOCATE(vec)
					ALLOCATE(vec(SIZE(vec_act)))
					vec = vec_act
				END IF
				
				! Re-normalizing VEC so that the numbers do not become too small
				vec = vec / DSQRT(SUM(vec**2))
				
				! Dellocating vec_act
				DEALLOCATE(vec_act)

			END DO

			! Updating accuracy
			accuracy  = SUM(ABS(img-bot_half))/392.0

		END SUBROUTINE GENERATE_FROM_BOT

END MODULE UTILITY