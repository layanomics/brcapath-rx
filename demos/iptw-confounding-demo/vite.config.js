import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Set base to './' for relative asset paths — works for both
// GitHub Pages project subdirectory and local file:// preview.
export default defineConfig({
  plugins: [react()],
  base: '/brcapath-rx/confounding-demo/',
})
