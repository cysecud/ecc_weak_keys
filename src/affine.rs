
#[derive(Clone, Copy, Debug,PartialEq,PartialOrd, Eq, Ord)]
pub struct AffinePoint<T>{
    pub x: T,
    pub y: T,

    /// Is this point the point at infinity? 0 = no, 1 = yes
    ///
    pub infinity: u8,
}
impl <T>AffinePoint<T>{
    /// Create a new AffinePoint with the given coordinates.
    pub fn new(x: T, y: T) -> Self {
        Self { x, y, infinity: 0 }
    }
    /* pub fn negate(&self) -> Self {
        AffinePoint {
            x: self.x,
            y: self.y.negate(),
            infinity: self.infinity,
        }
    } */
}

