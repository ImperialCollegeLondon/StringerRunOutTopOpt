!DEC$ FREEFORM

module mathUtils_mod
    use abaqusXIT_mod
    use iso_fortran_env
    implicit none 

    real(real64), parameter :: pi = 3.141592653589793
    real(real64), parameter, dimension(3,3) :: Imat = &
        reshape((/1.0d0, 0.0d0, 0.0d0, &
                  0.0d0, 1.0d0, 0.0d0, &
                  0.0d0, 0.0d0, 1.0d0 /), (/3,3/))

contains 


    function vecProd(a) result(p)
        !! Computes a product of the vector elements
        implicit none 

        ! output 
        real(real64) :: p 
        ! arguments 
        real(real64), intent(in) :: a(:)
        ! vars
        integer :: i

        ! vars init
        p = 1.0d0
        i = 0

        do i = 1, size(a)
            p = p * a(i)
        end do
    
    end function

    function vec4mul(a, b) result(p)
        !! computes the matrix product of two vectors
        implicit none 

        ! output
        real(real64) :: p(4,4)
        ! arguments
        real(real64), intent(in) :: a(4), b(4)
        ! vars
        integer :: i, j 

        ! vars init
        i = 0
        j = 0
        p = 0.0d0

        do i = 1, 4 
            do j = 1, 4
                p(i,j) = a(i) * b(j)
            end do
        end do
    end function


    function vec8mul(a, b) result(p)
        !! computes the matrix product of two vectors
        implicit none 

        ! output
        real(real64) :: p(8,8)
        ! arguments
        real(real64), intent(in) :: a(8), b(8)
        ! vars
        integer :: i, j 

        ! vars init
        i = 0
        j = 0
        p = 0.0d0

        do i = 1, 8 
            do j = 1, 8
                p(i,j) = a(i) * b(j)
            end do
        end do
    end function

    function vecNorm(a) result(norm)
        !! Computes the norm of a 3D vector
        implicit none

        ! output
        real(real64) :: norm 
        ! arguments
        real(real64), intent(in) :: a(:)

        ! vars init
        norm = 0.0d0

        norm = sqrt(sum(a**2))
    end function vecNorm

    function vec3cross(a, b) result(c)
        !! 3d vector cross product
        implicit none

        ! output
        real(real64) :: c(3)
        ! arguments
        real(real64) :: a(3), b(3)

        ! vars init
        c = 0.0d0

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end

    function mat2det(A) result(det)
        !! Performs a direct calculation of the determinant of a 2x2 matrix.
        implicit none 

        real(real64) :: det 
        real(real64), intent(in) :: A(2,2)

        ! vars init
        det = 0.0d0

        det = A(1,1)*A(2,2) - A(1,2) * A(2,1)
    end function

    function mat2inv(A) result(B)
        !! Performs a direct calculation of the inverse of a 2x2 matrix.
        implicit none 

        real(real64) :: B(2,2)
        real(real64), intent(in) :: A(2,2)

        real(real64) :: detinv 

        ! vars init
        B       = 0.0d0
        detinv  = 0.0d0

        detinv = 1.0d0 / mat2det(A)

        B(1,1) = detinv * A(2,2)
        B(2,2) = detinv * A(1,1)
        B(1,2) = -detinv * A(1,2)
        B(2,1) = -detinv * A(2,1)
    end function

    function mat3det(A) result(det)
        !! Performs a direct calculation of the inverse of a 3x3 matrix.
        implicit none

        ! output
        real(real64) :: det
        ! arguments
        real(real64), intent(in) :: A(3,3)

        ! vars init
        det = 0.0d0

        det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)   &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)   &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
    end function mat3det

    function mat3inv(A, det) result(B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        implicit none
        
        ! output
        real(real64)             :: B(3,3)
        ! arguments
        real(real64), intent(in) :: A(3,3)
        real(real64), intent(in) :: det
        ! vars
        real(real64)             :: detinv

        ! vars init
        B       = 0.0d0
        detinv  = 0.0d0

        detinv = 1.0d0 / det
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function mat3inv

    function matLump(M,n) result(A)
        !! lumps a matrix of size [n x n]
        implicit none 

        integer, intent(in) :: n 
        real(real64), intent(in) :: M(n, n)
        real(real64) :: A(n, n)

        integer :: i 

        ! vars init
        A = 0.0d0
        i = 0

        do i = 1, n
            A(i, i) = 1.0d0 * sum(M(i,:))
        end do
    end function

    function areaTriangle(a, b, c) result(area)
        !! computes the area of an arbitrary triangle
        implicit none 

        real(real64) :: area 
        real(real64), intent(in) :: a(2), b(2), c(2) 

        real(real64) :: mat(3,3), det

        ! vars init
        area    = 0.0d0
        mat     = 0.0d0
        det     = 0.0d0

        mat(1,:) = (/a(1), a(2), 1.0d0/)
        mat(2,:) = (/b(1), b(2), 1.0d0/)
        mat(3,:) = (/c(1), c(2), 1.0d0/)

        det = mat3det(mat)
        area = 0.5d0 * det

        if (area <= 0.0d0) then 
            call sroXIT(" >>>>>> mathUtils | areaTriangle <<<<<< negative area ")
        end if
    end function

    
    function areaQuadrilateral(a,b,c,d) result(area)
        !! computes the area of an arbitrary quadrilateral
        implicit none 

        real(real64) :: area 
        real(real64), intent(in) :: a(2), b(2), c(2), d(2)

        ! vars init
        area = 0.0d0

        area = area + areaTriangle(a,b,c)
        area = area + areaTriangle(a,c,d)
    end function

    function areaTriangle3D(a, b, c) result(area)
        !! computes the area of an arbitrary triangle
        implicit none 

        real(real64) :: area 
        real(real64), intent(in) :: a(3), b(3), c(3) 


        ! vars init
        area    = 0.0d0

        area = vecNorm(vec3cross(b - a, c - a)) * 0.5d0

        if (area <= 0.0d0) then 
            call sroXIT(" >>>>>> mathUtils | areaTriangle3D <<<<<< negative area ")
        end if
    end function

    
    function areaQuadrilateral3D(a,b,c,d) result(area)
        !! computes the area of an arbitrary quadrilateral
        implicit none 

        real(real64) :: area 
        real(real64), intent(in) :: a(3), b(3), c(3), d(3)

        ! vars init
        area = 0.0d0

        area = area + areaTriangle3D(a,b,c)
        area = area + areaTriangle3D(a,c,d)
    end function


    function volTetrahedron(a, b, c, d) result(vol)
        !! computes the volume of an arbitrary tetrahedron
        implicit none 

        ! output
        real(real64) :: vol
        ! arguments
        real(real64), intent(in) :: a(3), b(3), c(3), d(3)

        ! vars init
        vol = 0.0d0

        vol = (1.0d0 / 6.0d0) * abs(dot_product((a - d), vec3cross((b - d), (c - d))))
    end function volTetrahedron

    function volHexahedron(coords) result(vol)
        !! computes the volume of an arbitrary hexahedron
        implicit none 

        ! output
        real(real64) :: vol 
        ! arguments
        real(real64), intent(in) :: coords(8,3) 
        ! vars
        real(real64) :: v1(3), v2(3), v3(3), v4(3), v5(3), v6(3), v7(3), v8(3)

        ! vars init
        vol = 0.0d0
        v1  = 0.0d0
        v2  = 0.0d0
        v3  = 0.0d0
        v4  = 0.0d0
        v5  = 0.0d0
        v6  = 0.0d0
        v7  = 0.0d0
        v8  = 0.0d0

        v1 = coords(1, :)
        v2 = coords(2, :)
        v3 = coords(3, :)
        v4 = coords(4, :)
        v5 = coords(5, :)
        v6 = coords(6, :)
        v7 = coords(7, :)
        v8 = coords(8, :)

        vol = volTetrahedron(v1, v2, v3, v6) &
            + volTetrahedron(v1, v3, v4, v8) &
            + volTetrahedron(v6, v7, v8, v3) &
            + volTetrahedron(v6, v8, v5, v4) &
            + volTetrahedron(v1, v6, v3, v8)
    end function volHexahedron

    function volWedge(coords) result(vol)
        !! computes the volume of an arbitrary triangular prism
        implicit none

        real(real64) :: vol
        real(real64), intent(in) :: coords(6,3)

        real(real64) :: v1(3), v2(3), v3(3), v4(3), v5(3), v6(3)

        ! vars init
        vol = 0.0d0
        v1  = 0.0d0
        v2  = 0.0d0
        v3  = 0.0d0
        v4  = 0.0d0
        v5  = 0.0d0
        v6  = 0.0d0

        v1 = coords(1, :)
        v2 = coords(2, :)
        v3 = coords(3, :)
        v4 = coords(4, :)
        v5 = coords(5, :)
        v6 = coords(6, :)

        vol = volTetrahedron(v1, v5, v6, v4) + &
              volTetrahedron(v1, v2, v3, v5) + &
              volTetrahedron(v1, v3, v5, v6)
    end function

    function getPointCoords(coords, Nmat, n) result(c)
        implicit none 

        integer, intent(in) :: n 
        real(real64), intent(in) :: coords(:,:), Nmat(:,:)
        real(real64) :: c(n)

        integer :: i 

        ! vars init
        c = 0.0d0
        i = 0

        do i = 1, n 
            c(i) = dot_product(Nmat(1,:), coords(:,i))
        end do
    end function

    function vec3distance(a, b) result(d)
        implicit none 

        real(real64) :: d 
        real(real64), intent(in) :: a(:), b(:)
        integer :: i, n
        real(real64) :: aux

        ! vars init
        d   = 0.0d0
        i   = 0
        n   = 0
        aux = 0.0d0

        n = size(a)

        d = 0.0d0
        do i = 1, n 
            aux = a(i) - b(i)
            d = d + (aux * aux)
        end do
        d = sqrt(d)

    end function

end module