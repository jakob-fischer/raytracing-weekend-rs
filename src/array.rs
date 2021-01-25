    pub struct Array2d<T> {
        width: usize,
        height: usize,
        vec: Vec<T>,
    }
    
    impl<T: Copy> Array2d<T> {
        pub fn new(width: usize, height: usize, initial: &T) -> Self {
            Array2d::<T> {
                width,
                height,
                vec: (0..(width * height)).map(|_| *initial).collect(),
            }
        }
    }
    
    impl<T> Array2d<T> {
        pub fn get_width(&self) -> usize {
            self.width
        }
    
        pub fn get_height(&self) -> usize {
            self.height
        }
    
        pub fn get(&self, x: usize, y: usize) -> &T {
            self.vec.get(x * self.height + y).unwrap()
        }
    
        pub fn get_mut(&mut self, x: usize, y: usize) -> &mut T {
            self.vec.get_mut(x * self.height + y).unwrap()
        }
    }
