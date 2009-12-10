;; der compiler meckert wenn ich nicht alle variablen nutze. hoffentlich schmeisst er sie auch raus
(defmacro img-do ((img &key (centered nil) (centered-x nil) (centered-y nil)) &body body)
  `(let* ((w (array-dimension ,img 0))
	  (h (array-dimension ,img 1))
	  (w/2 (* w .5d0))
	  (h/2 (* h .5d0))
	  (/w (/ 1d0 w))
	  (/h (/ 1d0 h)))
    (declare (type (simple-array double-float (* *)) ,img))
    (declare (double-float w/2 h/2 /w /h))
    (dotimes (j h)
      (declare (fixnum j))
      (let ((y (* /h (if (or ,centered ,centered-y)
			     (- j h/2)
			     j))))
	(declare (double-float y))
	(dotimes (i w)
	  (declare (fixnum i))
	  (let ((x (* /w (if (or ,centered ,centered-x)
			     (- i w/2)
			     i))))
	    (declare (double-float x))
	    ,@body))))))

(defun extrema (img)
  (declare (type (simple-array double-float (* *)) img))
  (let* ((mi (aref img 0 0))
	 (ma mi))
    (declare (double-float mi ma))
    (img-do (img)
      (let ((v (aref img i j)))
	(setf mi (min mi v))
	(setf ma (max ma v))))
    (values mi ma)))

(defun normalize (img)
  (declare (type (simple-array double-float (* *)) img))
  (multiple-value-bind (mini maxi)
      (extrema img)
    (declare (double-float mini maxi))
    (let ((s (/ (- maxi mini))))
      (img-do (img)
	(setf (aref img i j)
	      (* s (- (aref img i j) mini))))))
  img)


(defun write-pgm (filename a &key (maxi 1d0 maxi-p) (mini 0d0 mini-p))
  "Write scalar real data from A into FILENAME as raw Portable Gray
Map image. Linearly scale values so that :MINI is 0 and :maxi is 255.
If those keyword arguments are not present the extrema of the data are
chosen."
  (declare (type (simple-array double-float (* *)) a))
  (let* ((h (array-dimension a 1))
         (w (array-dimension a 0))
         (buf (make-array (* w h) :element-type '(unsigned-byte 8))))
    (declare (type (simple-array (unsigned-byte 8) *) buf))
    (labels ((clamp (v) ;; cut value v so that it fits into [0..255]
	       (declare (fixnum v))
	       (cond
		 ((> v 255) 255)
		 ((< v 0) 0)
		 (t (the (unsigned-byte 8) v))))
	     (copy-to-buf (buf a mini maxi)
	       (declare (type (simple-array (unsigned-byte 8) *) buf))
	       (declare (type (simple-array double-float (* *)) a))
	       (declare (type double-float mini maxi))
	       (let* ((scale (* 255 (/ (- maxi mini)))))
		 (img-do (a)
		   (setf (aref buf (+ i (* w j)))
			 (clamp (floor (* scale (- (aref a i j) mini)))))))))
      (if (and maxi-p mini-p)
	  (copy-to-buf buf a mini maxi) ;; use scaling as provided by user
	  (multiple-value-bind (mini-b maxi-b) ;; find correct scaling from data
	      (extrema a)
	    (declare (double-float maxi-b mini-b))
	    (setf mini mini-b
		  maxi maxi-b)
	    (copy-to-buf buf a mini maxi))))
    
    (with-open-file (s filename
                       :direction :output
                       :if-exists :supersede
                       :element-type :default)
      (format s "P5 ~d ~d 255~%" w h)
      (write-sequence buf s)))
  (values maxi mini))


(defparameter lut (make-array '(256 256) :element-type 'double-float))

(img-do (lut)
  (setf (aref lut i j) (+ x y)))

(let* ((ll 16d0)
       (l .480))
  (labels ((uce->d (x) ; convert voltage U_CE=x V into mirror deflection d/um
	     (let* ((x3 (* x x x))
		    (a 0.0113608)
		    (d 4.31718))
	       (* 1/1000 (+ (* a x3) d)))))
    (let* ((imagefactor .75)
	   (maxvoltage 30)
	   (max-deflection (uce->d (* imagefactor maxvoltage)))
	   (max-piston-intensity (- 1 (cos (* 2 2 pi (/ l)
					      max-deflection)))))
        (format t "~s~%" (list "max-defl" max-deflection))  
      (img-do (lut)
	(let* ((d1 (uce->d (* x imagefactor maxvoltage)))
	       (d2 (uce->d (* y imagefactor maxvoltage)))
	       (delta-d (- d1 d2))
	       (arg (* 2 pi delta-d (/ l)))
	       (torsion-mirror-intensity (- 1 (if (= arg 0d0)
					   1d0
					   (/ (sin arg) arg)))))
	(setf (aref lut i j) (/ torsion-mirror-intensity max-piston-intensity)))))
    
    ))
  

(write-pgm "/dev/shm/lut.pgm" lut)


(defun read-without-comment (s)
  "Read the next data from stream S, skipping any comments."
  (read s)
;; (if (eq (peek-char nil s) #\#) ;; FIXME that should work
;;       (progn
;;         (read s)
;;         (read-without-comment s))
;;       (read s))
)

(defun read-pgm (filename)
  "Read binary PGM file."
  (with-open-file (s filename :external-format :ascii) 
    (when (not (equal (symbol-name (read s)) "P5"))
      (error "not a binary PGM file"))
    (let* ((width (read-without-comment s))
	   (height (read-without-comment s))
	   (grays (read-without-comment s))
	   (pos (file-position s))
	   (data (make-array (* width height) :element-type '(unsigned-byte 8)))
	   (bild (make-array (list width height) :element-type 'double-float)))
      (with-open-file (s filename :element-type '(unsigned-byte 8))
	(file-position s pos)
	(read-sequence data s))
      (img-do (bild)
	(setf (aref bild i j)  (coerce (aref data (+ i (* w j))) 'double-float)))
      bild)))


(defparameter erika (normalize (read-pgm "erika.pgm")))

; i want some really black region
(img-do (erika)
  (when (< (let ((x1 (- x .57))
		 (y1 (- y .35)))
	     (+ (* x1 x1) (* y1 y1))) .001)
    (setf (aref erika i j) 0d0)))

(defparameter erika2 (make-array '(256 256) :element-type 'double-float))

(defun get-lut (intensity dic-neighbor-val)
  "scan through the lookup table and look for the value that the
actual pixel should have so that the intensity goal is achieved
together with its predecessor in x direction."
  (labels ((merit (lut-pos)
	     (let* ((lut-intensity (aref lut lut-pos dic-neighbor-val))
		    (diff (- intensity lut-intensity)))
	       (* diff diff))))
    (let* ((lut-pos 0)
	   (mi (merit lut-pos))) 
      (dotimes (i 256) ;; FIXME: slowest possible algorithm :-(
	(let ((v (merit i)))
	  (when (< v mi)
	    (setf lut-pos i
		  mi v))))
      lut-pos)))


;; set leftmost pixel to zero angle then update neighbors on the right
;; (along x) according to the lookup table
(time 
 (img-do (erika2)
   (if (eq i 0)
       (setf (aref erika2 i j) 0d0)
       (setf (aref erika2 i j) 
	     (coerce (get-lut (aref erika i j)
			      (coerce (floor (aref erika2 (1- i) j)) 'fixnum))
		     'double-float)))))

(write-pgm "/dev/shm/erika.pgm" erika)
(write-pgm "/dev/shm/erika2.pgm" erika2)

(defparameter bmp-header
  #(#x42 #x4d #x36 #x00 #x03 #x00 #x00 #x00  #x00 #x00 #x36 #x00 #x00 #x00 #x28 #x00
    #x00 #x00 #x00 #x01 #x00 #x00 #x00 #x01  #x00 #x00 #x01 #x00 #x18 #x00 #x00 #x00
    #x00 #x00 #x00 #x00 #x03 #x00 #xc4 #x0e  #x00 #x00 #xc4 #x0e #x00 #x00 #x00 #x00
    #x00 #x00 #x00 #x00 #x00 #x00))

(let ((m 0))
  (img-do (erika2)
    (setf m (max m (aref erika2 i j))))
  m)


(defparameter erika3 (make-array (* 256 256 3) :element-type '(unsigned-byte 8)))

(img-do (erika2)
	(dotimes (k 3)
	  (setf (aref erika3 (+ k (* 3 (+ i (* w j)))))
		(coerce (floor (aref erika2 i j)) '(unsigned-byte 8)))))

(defun write-bmp24 (filename bgr)
  (with-open-file (s filename
		     :direction :output
		     :if-exists :supersede
		     :element-type :default)
    (write-sequence bmp-header s)
    (write-sequence bgr s)
    nil))

(write-bmp24 "/dev/shm/erika3.bmp" erika3)

(defparameter checkerboard (make-array '(256 256) :element-type 'double-float))


;(defmacro xor (v1 v2)
;`(not (eq (not ,v1) (not ,v2)))) 

;; this will generate exactly the same file as Schachbrett.bmp
(progn
  (dotimes (i (/ 256 8))
    (dotimes (j (/ 256 8))
      (dotimes (ii 8)
	(dotimes (jj 8)
	  (setf (aref checkerboard (+ (* 8 i) ii) (+ (* 8 j) jj))
		(if (= 0 (mod (+ i j) 2))
		    255d0
		    0d0))))))
  
  (defparameter cb3 (make-array (* 256 256 3) :element-type '(unsigned-byte 8)))
  
  (img-do (checkerboard)
    (dotimes (k 3)
      (setf (aref cb3 (+ k (* 3 (+ i (* w j)))))
	    (coerce (floor (aref checkerboard i j)) '(unsigned-byte 8)))))
  
  (write-bmp24 "/dev/shm/checkerboard3.bmp" cb3))


(defun apply-lut (image)
  (let ((image2 (make-array '(256 256) :element-type 'double-float)))
    (img-do (image2)
      (if (eq i 0)
	  (setf (aref image2 i j) 0d0)
	  (setf (aref image2 i j) 
		(coerce (get-lut (aref image i j)
				 (coerce (floor (aref image2 (1- i) j)) 'fixnum))
			'double-float)))
      nil)
    image2))
 

(defun gray2d->bgr1d (image)
  (let ((out (make-array (* 256 256 3) :element-type '(unsigned-byte 8))))
    (img-do (image)
      (dotimes (k 3)
	(setf (aref out (+ k (* 3 (+ i (* w j)))))
	      (coerce (floor (aref image i j)) '(unsigned-byte 8)))))))

(progn
  (let ((size 32) ;; must be power of two
	)
    (dotimes (i (/ 256 size))
    (dotimes (j (/ 256 size))
      (dotimes (ii size)
	(dotimes (jj size)
	  (setf (aref checkerboard (+ (* size i) ii) (+ (* size j) jj))
		(if (= 0 (mod (+ i j) 2))
		    (* 32d0 ii)
		    0d0)))))))
  
  
  (defparameter cb3 (make-array (* 256 256 3) :element-type '(unsigned-byte 8)))
  
  ;; (img-do (checkerboard)
;;     (dotimes (k 3)
;;       (setf (aref cb3 (+ k (* 3 (+ i (* w j)))))
;; 	    (coerce (floor (aref checkerboard i j)) '(unsigned-byte 8)))))
  
  ;;  (write-bmp24 "/dev/shm/gray-checkerboard.bmp" cb3)
  (write-pgm "/dev/shm/gray-checkerboard-with-lut.pgm" (apply-lut checkerboard))
  (write-bmp24 "/dev/shm/gray-checkerboard-with-lut.bmp" (gray2d->bgr1d (apply-lut checkerboard)))
  )

