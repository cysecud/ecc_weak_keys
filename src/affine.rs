use crate::field::FieldPoint;
#[derive(Clone, Copy, Debug,PartialEq,PartialOrd, Eq, Ord)]
pub struct AffinePoint {
    pub x: FieldPoint,
    pub y: FieldPoint,

    /// Is this point the point at infinity? 0 = no, 1 = yes
    ///
    pub infinity: u8,
}
impl AffinePoint {
    /// Create a new AffinePoint with the given coordinates.
    pub fn new(x: FieldPoint, y: FieldPoint) -> Self {
        Self { x, y, infinity: 0 }
    }
    pub fn negate(&self) -> Self {
        AffinePoint {
            x: self.x,
            y: self.y.negate(),
            infinity: self.infinity,
        }
    }
}

