#[derive(Debug, Clone,Eq,Hash, Copy,PartialEq)]
pub struct ProjectivePoint<T> {
    pub x: T,
    pub y: T,
    pub z: T,

    /// Is this point the point at infinity? 0 = no, 1 = yes
    ///
    pub infinity: u8,
}
impl<T> ProjectivePoint<T> {
    pub fn new(x: T, y: T,z:T) -> Self {
        Self { x, y, z,infinity: 0 }
    }  
    /* pub fn negate(&self) -> Self {
        ProjectivePoint  {
            x: self.x,
            y: self.y.negate(),
            z:self.z,
            infinity: self.infinity,
        }
    } */  
}
/* impl <T>Neg for ProjectivePoint<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        
    }
} */

